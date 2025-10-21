#!/usr/bin/env python3
"""
Refractory to Sensitive Shift Analysis
分析基因敲除是否能将耐药细胞状态推向敏感细胞状态。
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import pickle
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import multiprocessing as mp
from pathlib import Path
import scanpy as sc
from datasets import Dataset

# 修复导入：直接从geneformer库的顶层导入
from geneformer import EmbExtractor
from geneformer.in_silico_perturber import InSilicoPerturber
from geneformer.tokenizer import TranscriptomeTokenizer

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# =================================================================================
# == 自包含的辅助函数 (从Geneformer源码复制) ==
# =================================================================================

def load_model(model_type, num_classes, model_directory, mode="eval"):
    """
    直接加载模型，避免依赖perturber_utils中的复杂逻辑。
    这是一个简化的版本，专注于加载评估模型。
    """
    from transformers import BertForSequenceClassification, BertForMaskedLM

    model_classes = {
        "Pretrained": BertForMaskedLM,
        "GeneClassifier": BertForMaskedLM,  # 根据需要调整
        "CellClassifier": BertForSequenceClassification,
    }
    model_class = model_classes.get(model_type)
    if not model_class:
        raise ValueError(f"未知的模型类型: {model_type}")

    config_path = os.path.join(model_directory, 'config.json')
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"在 '{model_directory}' 中未找到 'config.json'。请确认模型路径是否正确。")

    model_args = {
        "pretrained_model_name_or_path": model_directory,
        "output_hidden_states": True,
        "output_attentions": False,
        "torch_dtype": torch.float16 if torch.cuda.is_available() else torch.float32,  # 混合精度
    }
    if model_type == "CellClassifier":
        model_args["num_labels"] = num_classes

    model = model_class.from_pretrained(**model_args)
    model.eval()
    
    # 性能优化
    if torch.cuda.is_available():
        model = model.to("cuda")
        model = model.half()  # 使用半精度以节省显存和提升速度
        # 启用CUDA优化
        torch.backends.cudnn.benchmark = True
    else:
        # CPU优化
        torch.set_num_threads(torch.get_num_threads())  # 使用所有可用CPU核心
    
    return model

def forward_pass_single_cell(model, example_cell, layer_to_quant):
    """
    对单个细胞进行前向传播。
    """
    # 获取input_ids并确保是正确的张量格式
    input_ids = example_cell["input_ids"]
    
    # 处理不同的数据格式
    if isinstance(input_ids, torch.Tensor):
        # 如果已经是张量，确保数据类型正确
        input_ids = input_ids.long()
    else:
        # 如果不是张量，转换为张量
        try:
            input_ids = torch.tensor(input_ids, dtype=torch.long)
        except:
            # 如果转换失败，尝试先转换为列表
            input_ids = torch.tensor(list(input_ids), dtype=torch.long)
    
    # 确保是2D张量 (batch_size, sequence_length)
    if input_ids.dim() == 1:
        input_ids = input_ids.unsqueeze(0)
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    
    with torch.no_grad():
        outputs = model(input_ids=input_ids.to(device))
    
    # 确保hidden_states在预期的位置
    if isinstance(outputs, dict):
        hidden_states = outputs.get('hidden_states')
    else:
        # 兼容元组输出格式
        hidden_states = outputs[1]
        
    if hidden_states is None:
        raise ValueError("模型输出中未找到 'hidden_states'。")
        
    emb = torch.squeeze(hidden_states[layer_to_quant])
    del outputs
    return emb

def mean_nonpadding_embs(embs, original_lens, dim=1):
    """
    计算非填充部分的平均嵌入。
    """
    # 确保embs是3D张量
    if embs.dim() == 2:
        embs = embs.unsqueeze(0)
        
    # 确保original_lens是张量
    if not isinstance(original_lens, torch.Tensor):
        original_lens = torch.tensor(original_lens, device=embs.device)
        
    mask = torch.arange(embs.size(dim), device=embs.device) < original_lens.unsqueeze(1)
    
    masked_embs = embs.masked_fill(~mask.unsqueeze(2), 0.0)
    mean_embs = masked_embs.sum(dim) / original_lens.view(-1, 1).float()
    
    return mean_embs

# =================================================================================
# == 主要分析脚本 ==
# =================================================================================

def setup_logging(results_dir):
    """设置日志记录，同时打印到控制台和文件"""
    log_file_path = os.path.join(results_dir, "analysis_log.txt")
    
    def log_message(message, also_print=True):
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        formatted_message = f"[{timestamp}] {message}"
        if also_print:
            print(formatted_message)
        with open(log_file_path, 'a', encoding='utf-8') as f:
            f.write(formatted_message + '\n')
    
    return log_message

def load_gene_mapping(log_func, tokenizer):
    """加载基因名称映射字典"""
    log_func("开始加载基因名称映射字典...")
    project_root = Path(__file__).resolve().parent.parent
    gene_name_dict_path = project_root / "geneformer" / "gene_name_id_dict_gc104M.pkl"
    
    try:
        with open(gene_name_dict_path, 'rb') as f:
            gene_name_id_dict = pickle.load(f)
        
        id_to_name_dict = {v: k for k, v in gene_name_id_dict.items()}
        name_to_id_dict = gene_name_id_dict
        
        log_func(f"成功加载基因名称字典，包含 {len(id_to_name_dict)} 个基因。")
        return id_to_name_dict, name_to_id_dict
    except Exception as e:
        log_func(f"错误：加载基因名称字典失败 - {e}")
        return {}, {}

def get_embeddings_batch(model, dataset, batch_size=32, layer_to_quant=-1):
    """批量提取嵌入向量以提升性能"""
    embeddings = []
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    # 添加进度条
    total_batches = (len(dataset) + batch_size - 1) // batch_size
    
    for i in tqdm(range(0, len(dataset), batch_size), 
                  desc=f"批量处理嵌入提取", 
                  total=total_batches):
        batch_end = min(i + batch_size, len(dataset))
        batch_indices = list(range(i, batch_end))
        batch_dataset = dataset.select(batch_indices)
        
        # 准备批量输入
        batch_input_ids = []
        batch_lengths = []
        
        for j in range(len(batch_dataset)):
            input_ids = batch_dataset[j]["input_ids"]
            if isinstance(input_ids, torch.Tensor):
                input_ids = input_ids.long()
            else:
                try:
                    input_ids = torch.tensor(input_ids, dtype=torch.long)
                except:
                    input_ids = torch.tensor(list(input_ids), dtype=torch.long)
            
            batch_input_ids.append(input_ids)
            batch_lengths.append(len(input_ids))
        
        # 填充到相同长度并创建attention_mask
        max_len = max(batch_lengths)
        padded_inputs = []
        attention_masks = []
        
        for input_ids, length in zip(batch_input_ids, batch_lengths):
            if len(input_ids) < max_len:
                padding = torch.zeros(max_len - len(input_ids), dtype=torch.long)
                input_ids = torch.cat([input_ids, padding])
            padded_inputs.append(input_ids)
            
            # 创建attention mask (1表示真实token，0表示padding)
            mask = torch.ones(max_len, dtype=torch.long)
            mask[length:] = 0
            attention_masks.append(mask)
        
        # 转换为批量张量
        batch_tensor = torch.stack(padded_inputs).to(device)
        attention_mask = torch.stack(attention_masks).to(device)
        
        with torch.no_grad():
            outputs = model(input_ids=batch_tensor, attention_mask=attention_mask)
            if isinstance(outputs, dict):
                hidden_states = outputs.get('hidden_states')
            else:
                hidden_states = outputs[1]
            
            if hidden_states is None:
                raise ValueError("模型输出中未找到 'hidden_states'。")
            
            # 提取最后一层的嵌入
            batch_embs = hidden_states[layer_to_quant]  # [batch_size, seq_len, hidden_dim]
            
            # 对每个样本计算平均嵌入（忽略填充部分）
            for j, length in enumerate(batch_lengths):
                emb = batch_embs[j, :length, :].mean(dim=0)  # 平均非填充token
                embeddings.append(emb)
        
        # 清理GPU内存
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    
    return torch.stack(embeddings)

def get_sensitive_avg_embedding(model, sensitive_dataset, emb_dir, log_func):
    """计算或加载敏感组的平均嵌入（优化版本）"""
    avg_emb_path = os.path.join(emb_dir, "sensitive_avg_embedding.pt")
    
    if os.path.exists(avg_emb_path):
        log_func(f"从缓存加载敏感组平均嵌入: {avg_emb_path}")
        return torch.load(avg_emb_path)

    log_func(f"计算 {len(sensitive_dataset)} 个敏感组细胞的嵌入...")
    
    # 使用批处理提取嵌入（性能优化）
    batch_size = 8 if torch.cuda.is_available() else 4  # 调整批大小避免内存溢出
    log_func(f"使用批大小: {batch_size}")
    
    all_embs = get_embeddings_batch(model, sensitive_dataset, batch_size, layer_to_quant=-1)
    avg_embedding = torch.mean(all_embs, dim=0)
    
    log_func(f"敏感组平均嵌入计算完成，形状: {avg_embedding.shape}")
    
    torch.save(avg_embedding, avg_emb_path)
    log_func(f"敏感组平均嵌入已缓存: {avg_emb_path}")
    
    return avg_embedding

def get_knockout_embeddings(model, cell_dataset, gene_token_to_knockout):
    """为单个细胞生成敲除前后的嵌入"""
    
    original_emb_full = forward_pass_single_cell(model, cell_dataset, layer_to_quant=-1)
    original_cell_emb = mean_nonpadding_embs(original_emb_full, cell_dataset["length"])

    input_ids = cell_dataset["input_ids"][0]
    knockout_ids = [token for token in input_ids if token != gene_token_to_knockout]
    
    knockout_dataset = Dataset.from_dict({
        "input_ids": [knockout_ids],
        "length": [len(knockout_ids)]
    })
    
    knockout_emb_full = forward_pass_single_cell(model, knockout_dataset, layer_to_quant=-1)
    knockout_cell_emb = mean_nonpadding_embs(knockout_emb_full, knockout_dataset["length"])

    return original_cell_emb.squeeze(), knockout_cell_emb.squeeze()

def tokenize_anndata(adata_path, output_dir, log_func):
    """对AnnData文件进行分词"""
    tokenized_dataset_path = output_dir / "tokenized_data.dataset"
    if tokenized_dataset_path.exists():
        log_func(f"已找到分词后的数据集，从缓存加载: {tokenized_dataset_path}")
        from datasets import load_from_disk
        return load_from_disk(str(tokenized_dataset_path))

    log_func(f"开始对 {adata_path} 进行分词...")
    from geneformer.tokenizer import TranscriptomeTokenizer
    
    # 确保输出目录存在
    output_dir.mkdir(exist_ok=True)
    
    tokenizer = TranscriptomeTokenizer(nproc=8)
    tokenizer.tokenize_data(
        data_directory=os.path.dirname(adata_path),
        output_directory=str(output_dir),
        output_prefix="tokenized_data",
        file_format="h5ad"
    )
    
    # 加载刚刚创建的数据集
    from datasets import load_from_disk
    tokenized_dataset = load_from_disk(os.path.join(output_dir, "tokenized_data.dataset"))
    
    log_func(f"分词完成！数据已保存到: {output_dir}")
    return tokenized_dataset

def main():
    # 性能优化配置
    torch.set_num_threads(min(mp.cpu_count(), 16))  # 限制CPU线程数避免过度占用
    if torch.cuda.is_available():
        torch.backends.cudnn.benchmark = True  # 优化CUDA性能
        torch.backends.cudnn.deterministic = False  # 启用非确定性算法以提升速度
        print(f"检测到GPU: {torch.cuda.get_device_name()}")
        print(f"GPU内存: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")
    else:
        print(f"使用CPU，线程数: {torch.get_num_threads()}")
    
    # 1. 设置路径和日志
    script_name = Path(__file__).stem
    project_root = Path(__file__).resolve().parent.parent
    results_dir = project_root / 'results' / script_name
    results_dir.mkdir(exist_ok=True)
    
    emb_dir = results_dir / 'embeddings'
    emb_dir.mkdir(exist_ok=True)
    
    log_func = setup_logging(results_dir)
    log_func("="*60)
    log_func("耐药细胞向敏感细胞状态转变分析开始")
    log_func(f"所有输出将保存在: {results_dir}")
    log_func("="*60)

    # 2. 加载和准备数据
    try:
        log_func("加载/分词数据...")
        project_root = Path(__file__).resolve().parent.parent
        adata_path = project_root / "data" / "processed" / "geneformer_input.h5ad"
        tokenized_data_dir = results_dir / "tokenized"
        
        # 步骤1: 分词
        tokenized_dataset = tokenize_anndata(adata_path, tokenized_data_dir, log_func)
        
        # 步骤2: 加载模型和映射
        log_func("加载模型和基因映射...")
        model_path = project_root / "Geneformer-V2-104M_CLcancer"
        tokenizer = TranscriptomeTokenizer() # 再次创建以获取字典
        
        model = load_model("CellClassifier", 2, str(model_path))
        id_to_name, name_to_id = load_gene_mapping(log_func, tokenizer)

        # 步骤3: 拆分数据集
        # 使用原始adata中的status列来过滤tokenized_dataset
        original_adata = sc.read_h5ad(adata_path)
        sensitive_indices = [i for i, ct in enumerate(original_adata.obs['status']) if ct == 'Sensitive']
        refractory_indices = [i for i, ct in enumerate(original_adata.obs['status']) if ct == 'Refractory']

        sensitive_dataset = tokenized_dataset.select(sensitive_indices)
        refractory_dataset = tokenized_dataset.select(refractory_indices)
        
        log_func(f"数据集拆分完成: {len(sensitive_dataset)}个敏感细胞, {len(refractory_dataset)}个耐药细胞")

    except Exception as e:
        log_func(f"错误：加载模型或数据失败: {e}")
        import traceback
        log_func(traceback.format_exc())
        return

    # 3. 计算敏感组平均嵌入
    sensitive_avg_emb = get_sensitive_avg_embedding(model, sensitive_dataset, emb_dir, log_func)

    # 4. 遍历基因，进行分析
    genes_to_analyze = ["SMARCA4", "TP53", "PTEN", "BRCA1", "BRCA2", "SMARCA2", "ABCB1"]
    all_shift_data = []

    for gene_name in genes_to_analyze:
        gene_id = name_to_id.get(gene_name)
        if not gene_id:
            log_func(f"警告：在名称字典中未找到基因 {gene_name}")
            continue
            
        gene_token = tokenizer.gene_token_dict.get(gene_id)
        if not gene_token:
            log_func(f"警告：在token字典中未找到基因 {gene_id} ({gene_name})")
            continue

        log_func(f"\n--- 分析基因: {gene_name} (ID: {gene_id}, Token: {gene_token}) ---")

        # 从已经分词的Dataset中筛选
        target_refractory_dataset = refractory_dataset.filter(
            lambda example: gene_token in example['input_ids'], num_proc=4
        )
        
        if len(target_refractory_dataset) == 0:
            log_func("未找到表达该基因的耐药细胞。")
            continue
            
        log_func(f"找到 {len(target_refractory_dataset)} 个表达该基因的耐药细胞。开始逐个处理...")

        for i in tqdm(range(len(target_refractory_dataset)), desc=f"处理 {gene_name} 细胞"):
            cell_dataset = target_refractory_dataset.select([i])
            
            try:
                pre_ko_emb, post_ko_emb = get_knockout_embeddings(model, cell_dataset, gene_token)
                
                cos = torch.nn.CosineSimilarity(dim=0)
                pre_sim = cos(pre_ko_emb.cpu(), sensitive_avg_emb.cpu()).item()
                post_sim = cos(post_ko_emb.cpu(), sensitive_avg_emb.cpu()).item()
                
                shift = post_sim - pre_sim
                
                # 假设cell_id在原始adata的obs中
                original_cell_index = refractory_indices[i]
                cell_id = original_adata.obs_names[original_cell_index]

                all_shift_data.append({
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "cell_id": cell_id,
                    "pre_knockout_similarity": pre_sim,
                    "post_knockout_similarity": post_sim,
                    "similarity_shift": shift
                })
            except Exception as e:
                original_cell_index = refractory_indices[i]
                cell_id = original_adata.obs_names[original_cell_index]
                log_func(f"处理细胞 {cell_id} 时出错: {e}")

    # 5. 保存结果
    log_func("\n--- 所有基因处理完毕，正在保存结果... ---")
    shift_df = pd.DataFrame(all_shift_data)
    results_path = results_dir / "similarity_shift_results.csv"
    shift_df.to_csv(results_path, index=False)
    log_func(f"所有基因的相似度变化数据已保存: {results_path}")
    
    log_func("="*60)
    log_func("数据生成和计算阶段完成！下一步将进行可视化。")
    log_func("="*60)

if __name__ == "__main__":
    main()
