#!/usr/bin/env python3
"""
卵巢癌化疗耐药相关基因敲除分析
专注于7个关键基因：SMARCA2, SMARCA4, BRCA1, BRCA2, TP53, PTEN, MDR1
"""

import os
import sys
import gc
import torch
import pickle
import numpy as np
import pandas as pd
from datetime import datetime
import logging
from typing import Dict, List, Optional, Tuple
from datasets import Dataset, concatenate_datasets

# 设置环境变量
os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'

# 添加geneformer路径
sys.path.append('/f%3A/githubclone/Geneformer')

from geneformer import EmbExtractor, InSilicoPerturber
from geneformer.in_silico_perturber_stats import InSilicoPerturberStats

def setup_logging(output_dir: str) -> logging.Logger:
    """设置日志"""
    log_file = os.path.join(output_dir, f"gene_knockout_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger(__name__)

def clear_gpu_memory():
    """清理GPU内存"""
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()
    gc.collect()

class OvarianCancerGeneKnockoutAnalyzer:
    """卵巢癌基因敲除分析器"""
    
    def __init__(self, 
                 model_path: str = "Geneformer-V2-104M_CLcancer/",
                 output_dir: str = "results/targeted_gene_knockout",
                 max_ncells: int = 400,
                 emb_batch_size: int = 8,
                 perturb_batch_size: int = 4,
                 nproc: int = 2):
        
        self.model_path = model_path
        self.output_dir = output_dir
        self.max_ncells = max_ncells
        self.emb_batch_size = emb_batch_size
        self.perturb_batch_size = perturb_batch_size
        self.nproc = nproc
        
        # 目标基因列表 - 7个与卵巢癌化疗耐药相关的关键基因
        self.target_genes_symbols = [
            "SMARCA2",    # SWI/SNF相关基因，染色质重塑
            "SMARCA4",    # SWI/SNF相关基因，染色质重塑  
            "BRCA1",      # DNA修复，化疗敏感性
            "BRCA2",      # DNA修复，铂类耐药
            "TP53",       # 肿瘤抑制基因，细胞周期/凋亡
            "PTEN",       # PI3K/AKT通路，顺铂敏感性
            "ABCB1"       # 多药耐药基因1 (也称为MDR1)
        ]
        
        # 基因别名映射
        self.gene_aliases = {
            "MDR1": ["ABCB1", "PGY1", "GP170"],
            "SMARCA2": ["BRM", "SNF2L2"],
            "SMARCA4": ["BRG1", "SNF2L4"],
            "BRCA1": ["BRCAI", "BRCC1"],
            "BRCA2": ["BRCC2", "FANCD1"],
            "TP53": ["P53", "TRP53"],
            "PTEN": ["MMAC1", "TEP1"]
        }
        
        os.makedirs(output_dir, exist_ok=True)
        self.logger = setup_logging(output_dir)
        
        # 加载Geneformer的基因名称到ID字典
        self.load_gene_name_id_dict()
        
        # 将基因符号转换为Ensembl ID
        self.target_genes = []
        for symbol in self.target_genes_symbols:
            if symbol in self.gene_name_id_dict:
                ensembl_id = self.gene_name_id_dict[symbol]
                self.target_genes.append(ensembl_id)
                self.logger.info(f"基因 {symbol} -> {ensembl_id}")
            else:
                self.logger.warning(f"基因 {symbol} 未在基因名称字典中找到")
        
        if not self.target_genes:
            raise ValueError("没有找到任何有效的目标基因")
        
        self.logger.info("=== 卵巢癌化疗耐药基因敲除分析 ===")
        self.logger.info(f"目标基因: {', '.join([f'{symbol}({ensembl})' for symbol, ensembl in zip(self.target_genes_symbols, self.target_genes)])}")
        self.logger.info(f"最大细胞数: {max_ncells}")
        self.logger.info(f"嵌入批次大小: {emb_batch_size}")
        self.logger.info(f"扰动批次大小: {perturb_batch_size}")
    
    def load_gene_name_id_dict(self):
        """加载Geneformer的基因名称到ID字典"""
        try:
            # 导入Geneformer的基因名称字典文件路径
            from geneformer import ENSEMBL_DICTIONARY_FILE
            
            import pickle
            with open(ENSEMBL_DICTIONARY_FILE, 'rb') as f:
                self.gene_name_id_dict = pickle.load(f)
                
            self.logger.info(f"成功加载基因名称字典，包含 {len(self.gene_name_id_dict)} 个基因")
            
        except Exception as e:
            self.logger.error(f"加载基因名称字典失败: {e}")
            # 如果加载失败，使用备用映射
            self.gene_name_id_dict = {
                'SMARCA2': 'ENSG00000080503',
                'SMARCA4': 'ENSG00000127616', 
                'BRCA1': 'ENSG00000012048',
                'BRCA2': 'ENSG00000139618',
                'TP53': 'ENSG00000141510',
                'PTEN': 'ENSG00000171862',
                'ABCB1': 'ENSG00000085563'
            }
            self.logger.warning("使用备用基因名称映射")
        
    def load_and_prepare_data(self, refractory_path: str, sensitive_path: str) -> Tuple[Dataset, Dataset]:
        """加载和准备数据"""
        self.logger.info("加载数据...")
        
        try:
            # 加载数据集
            refractory_dataset = Dataset.load_from_disk(refractory_path)
            sensitive_dataset = Dataset.load_from_disk(sensitive_path)
            
            self.logger.info(f"耐药组细胞数: {len(refractory_dataset)}")
            self.logger.info(f"敏感组细胞数: {len(sensitive_dataset)}")
            
            # 限制细胞数量以适应内存限制
            if len(refractory_dataset) > self.max_ncells:
                refractory_dataset = refractory_dataset.select(range(self.max_ncells))
                self.logger.info(f"耐药组限制为: {len(refractory_dataset)} 个细胞")
            
            if len(sensitive_dataset) > self.max_ncells:
                sensitive_dataset = sensitive_dataset.select(range(self.max_ncells))
                self.logger.info(f"敏感组限制为: {len(sensitive_dataset)} 个细胞")
            
            return refractory_dataset, sensitive_dataset
            
        except Exception as e:
            self.logger.error(f"数据加载失败: {e}")
            raise
    
    def find_gene_in_dataset(self, dataset: Dataset, gene_name: str) -> Optional[int]:
        """在数据集中查找基因ID"""
        # 检查数据集中的基因映射
        if hasattr(dataset, 'features') and 'input_ids' in dataset.features:
            # 尝试从第一个样本获取基因信息
            sample = dataset[0]
            input_ids = sample['input_ids']
            
            # 这里需要根据Geneformer的具体实现来查找基因ID
            # 通常需要使用tokenizer的基因字典
            self.logger.info(f"正在查找基因 {gene_name} 的ID...")
            
        return None
    
    def extract_embeddings(self, dataset: Dataset, group_name: str) -> str:
        """提取细胞嵌入"""
        self.logger.info(f"提取 {group_name} 组的细胞嵌入...")
        
        clear_gpu_memory()
        
        try:
            # 配置嵌入提取器
            embex = EmbExtractor(
                model_type="CellClassifier",
                num_classes=2,
                max_ncells=len(dataset),
                forward_batch_size=self.emb_batch_size,
                nproc=self.nproc,
                emb_mode="cls",
                emb_layer=-1
            )
            
            # 先保存dataset到临时文件
            temp_dataset_path = os.path.join(self.output_dir, f"temp_{group_name}_dataset")
            dataset.save_to_disk(temp_dataset_path)
            
            # 提取嵌入
            output_path = os.path.join(self.output_dir, f"{group_name}_embeddings")
            os.makedirs(output_path, exist_ok=True)
            
            embs = embex.extract_embs(
                model_directory=self.model_path,
                input_data_file=temp_dataset_path,
                output_directory=output_path,
                output_prefix=f"{group_name}_emb"
            )
            
            self.logger.info(f"{group_name} 组嵌入提取完成")
            clear_gpu_memory()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"{group_name} 组嵌入提取失败: {e}")
            clear_gpu_memory()
            raise
    
    def perform_gene_knockout(self, dataset: Dataset, gene_name: str, group_name: str) -> str:
        """执行基因敲除"""
        self.logger.info(f"对 {group_name} 组执行 {gene_name} 基因敲除...")
        
        clear_gpu_memory()
        
        try:
            # 配置扰动器
            isp = InSilicoPerturber(
                perturb_type="delete",
                perturb_rank_shift=None,
                genes_to_perturb=[gene_name],
                model_type="CellClassifier", 
                num_classes=2,
                max_ncells=len(dataset),
                forward_batch_size=self.perturb_batch_size,
                nproc=self.nproc,
                emb_mode="cls",
                emb_layer=-1
            )
            
            # 先保存dataset到临时文件
            temp_dataset_path = os.path.join(self.output_dir, f"temp_{group_name}_{gene_name}_dataset")
            dataset.save_to_disk(temp_dataset_path)
            
            # 执行扰动
            output_path = os.path.join(self.output_dir, f"{group_name}_{gene_name}_knockout")
            os.makedirs(output_path, exist_ok=True)
            
            isp.perturb_data(
                model_directory=self.model_path,
                input_data_file=temp_dataset_path,
                output_directory=output_path,
                output_prefix=f"{group_name}_{gene_name}_ko"
            )
            
            self.logger.info(f"{group_name} 组 {gene_name} 基因敲除完成")
            clear_gpu_memory()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"{group_name} 组 {gene_name} 基因敲除失败: {e}")
            clear_gpu_memory()
            raise
    
    def analyze_knockout_effects(self, knockout_results: Dict) -> pd.DataFrame:
        """分析敲除效果"""
        self.logger.info("分析基因敲除效果...")
        
        results_summary = []
        
        # 遍历所有基因符号，而不是Ensembl ID
        for i, gene_symbol in enumerate(self.target_genes_symbols):
            for group in ['refractory', 'sensitive']:
                try:
                    result_key = f"{group}_{gene_symbol}"
                    result_path = knockout_results.get(result_key)
                    
                    if result_path and os.path.exists(result_path):
                        # 这里添加具体的分析逻辑
                        # 例如：计算表达变化、细胞状态变化等
                        
                        result_info = {
                            'gene': gene_symbol,
                            'ensembl_id': self.target_genes[i],
                            'group': group,
                            'status': 'completed',
                            'output_path': result_path
                        }
                        results_summary.append(result_info)
                        self.logger.info(f"✓ {group} 组 {gene_symbol} 基因敲除结果已记录")
                        
                    else:
                        result_info = {
                            'gene': gene_symbol,
                            'ensembl_id': self.target_genes[i],
                            'group': group,
                            'status': 'failed',
                            'error': f"结果文件不存在: {result_path}"
                        }
                        results_summary.append(result_info)
                        self.logger.warning(f"✗ {group} 组 {gene_symbol} 基因敲除结果文件缺失")
                        
                except Exception as e:
                    self.logger.error(f"分析 {group} 组 {gene_symbol} 基因敲除结果失败: {e}")
                    
                    result_info = {
                        'gene': gene_symbol,
                        'ensembl_id': self.target_genes[i] if i < len(self.target_genes) else 'Unknown',
                        'group': group,
                        'status': 'failed',
                        'error': str(e)
                    }
                    results_summary.append(result_info)
        
        results_df = pd.DataFrame(results_summary)
        
        # 保存结果摘要
        summary_file = os.path.join(self.output_dir, "knockout_results_summary.csv")
        results_df.to_csv(summary_file, index=False, encoding='utf-8')
        self.logger.info(f"结果摘要已保存到: {summary_file}")
        
        return results_df
    
    def run_analysis(self, refractory_path: str, sensitive_path: str):
        """运行完整分析"""
        self.logger.info("开始卵巢癌化疗耐药基因敲除分析...")
        
        try:
            # 1. 加载数据
            refractory_dataset, sensitive_dataset = self.load_and_prepare_data(
                refractory_path, sensitive_path
            )
            
            # 2. 提取基线嵌入
            self.logger.info("提取基线细胞嵌入...")
            refractory_emb_path = self.extract_embeddings(refractory_dataset, "refractory")
            sensitive_emb_path = self.extract_embeddings(sensitive_dataset, "sensitive")
            
            # 3. 执行基因敲除实验
            knockout_results = {}
            
            for i, gene_ensembl in enumerate(self.target_genes):
                gene_symbol = self.target_genes_symbols[i]
                self.logger.info(f"=== 开始 {gene_symbol} 基因敲除实验 ===")
                
                try:
                    # 耐药组敲除
                    refractory_ko_path = self.perform_gene_knockout(
                        refractory_dataset, gene_ensembl, "refractory"
                    )
                    knockout_results[f"refractory_{gene_symbol}"] = refractory_ko_path
                    
                    # 敏感组敲除
                    sensitive_ko_path = self.perform_gene_knockout(
                        sensitive_dataset, gene_ensembl, "sensitive"
                    )
                    knockout_results[f"sensitive_{gene_symbol}"] = sensitive_ko_path
                    
                    self.logger.info(f"{gene_symbol} 基因敲除实验完成")
                    
                except Exception as e:
                    self.logger.error(f"{gene_symbol} 基因敲除实验失败: {e}")
                    continue
            
            # 4. 分析结果
            results_summary = self.analyze_knockout_effects(knockout_results)
            
            # 5. 生成报告
            self.generate_report(results_summary)
            
            self.logger.info("=== 分析完成 ===")
            self.logger.info(f"结果保存在: {self.output_dir}")
            
        except Exception as e:
            self.logger.error(f"分析失败: {e}")
            raise
    
    def generate_report(self, results_summary: pd.DataFrame):
        """生成分析报告"""
        self.logger.info("生成分析报告...")
        
        report_file = os.path.join(self.output_dir, "analysis_report.md")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("# 卵巢癌化疗耐药基因敲除分析报告\n\n")
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## 分析概述\n\n")
            f.write(f"- 目标基因数量: {len(self.target_genes)}\n")
            f.write(f"- 目标基因: {', '.join(self.target_genes)}\n")
            f.write(f"- 最大细胞数: {self.max_ncells}\n\n")
            
            f.write("## 目标基因功能说明\n\n")
            gene_functions = {
                "SMARCA2": "SWI/SNF染色质重塑复合体成分，调节基因转录",
                "SMARCA4": "SWI/SNF染色质重塑复合体成分，调节基因转录", 
                "BRCA1": "DNA损伤修复，与化疗敏感性密切相关",
                "BRCA2": "DNA损伤修复，参与同源重组修复",
                "TP53": "肿瘤抑制基因，调控细胞周期和凋亡",
                "PTEN": "磷酸酶，调控PI3K/AKT信号通路",
                "MDR1": "多药耐药基因，编码P-糖蛋白转运蛋白"
            }
            
            for gene, function in gene_functions.items():
                f.write(f"- **{gene}**: {function}\n")
            
            f.write("\n## 实验结果摘要\n\n")
            
            # 统计成功和失败的实验
            completed = results_summary[results_summary['status'] == 'completed']
            failed = results_summary[results_summary['status'] == 'failed']
            
            f.write(f"- 成功完成的实验: {len(completed)}/{len(results_summary)}\n")
            f.write(f"- 失败的实验: {len(failed)}/{len(results_summary)}\n\n")
            
            if len(completed) > 0:
                f.write("### 成功完成的基因敲除实验\n\n")
                for _, row in completed.iterrows():
                    f.write(f"- {row['group']}组 {row['gene']} 基因敲除\n")
            
            if len(failed) > 0:
                f.write("\n### 失败的实验\n\n")
                for _, row in failed.iterrows():
                    f.write(f"- {row['group']}组 {row['gene']} 基因敲除: {row.get('error', '未知错误')}\n")
            
            f.write(f"\n## 输出文件位置\n\n")
            f.write(f"所有结果文件保存在: `{self.output_dir}`\n\n")
            
        self.logger.info(f"分析报告已保存到: {report_file}")

def main():
    """主函数"""
    # 设置GPU内存限制
    if torch.cuda.is_available():
        torch.cuda.set_per_process_memory_fraction(0.7)
        print(f"CUDA设备: {torch.cuda.get_device_name()}")
        print(f"可用显存: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")
    
    # 数据路径
    refractory_path = "data/tokenized_data/refractory_cells.dataset"
    sensitive_path = "data/tokenized_data/sensitive_cells.dataset"
    
    # 检查数据文件是否存在
    if not os.path.exists(refractory_path):
        print(f"错误: 找不到耐药组数据文件 {refractory_path}")
        return
    
    if not os.path.exists(sensitive_path):
        print(f"错误: 找不到敏感组数据文件 {sensitive_path}")
        return
    
    # 创建分析器并运行分析
    analyzer = OvarianCancerGeneKnockoutAnalyzer(
        model_path="Geneformer-V2-104M_CLcancer/",
        output_dir="results/targeted_gene_knockout",
        max_ncells=400,  # 减少细胞数以适应6GB显存
        emb_batch_size=8,
        perturb_batch_size=4,
        nproc=2
    )
    
    analyzer.run_analysis(refractory_path, sensitive_path)

if __name__ == "__main__":
    main()
