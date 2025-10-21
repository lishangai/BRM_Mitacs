import os
# 设置环境变量以避免 OpenMP 错误 - 必须在其他导入之前
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import pandas as pd
import scanpy as sc
import gseapy
import logging
import pickle

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('results/pathway_analysis.log'),
        logging.StreamHandler()
    ]
)

def load_gene_mapping_dict():
    """加载基因名到ID的映射字典"""
    dict_path = 'geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl'
    try:
        with open(dict_path, 'rb') as f:
            gene_dict = pickle.load(f)
        logging.info(f"成功加载基因映射字典，包含 {len(gene_dict)} 个条目")
        return gene_dict
    except Exception as e:
        logging.error(f"加载基因映射字典失败: {e}")
        return None

def create_ensembl_to_symbol_mapping(gene_dict):
    """创建从Ensembl ID到基因符号的映射"""
    if gene_dict is None:
        return None
    
    # 检查字典结构
    sample_items = list(gene_dict.items())[:5]
    logging.info(f"字典样本条目: {sample_items}")
    
    # 创建反向映射：从Ensembl ID到基因符号
    ensembl_to_symbol = {}
    for symbol, ensembl_id in gene_dict.items():
        if isinstance(ensembl_id, str) and ensembl_id.startswith('ENSG'):
            ensembl_to_symbol[ensembl_id] = symbol
    
    logging.info(f"创建了 {len(ensembl_to_symbol)} 个Ensembl ID到基因符号的映射")
    return ensembl_to_symbol

def main():
    logging.info("通路富集分析 (GSEA) 脚本开始运行。")
    
    # 定义目标基因和文件路径
    target_gene = 'SMARCA2'
    logging.info(f"将根据 '{target_gene}' 的敲除效果进行分组。")
    
    similarity_file = 'results/refractory_to_sensitive_shift_analysis/similarity_shift_results.csv'
    expression_file = 'data/processed/geneformer_input.h5ad'
    
    logging.info(f"相似度数据路径: {similarity_file}")
    logging.info(f"基因表达数据路径: {expression_file}")
    
    # 加载基因映射字典
    logging.info("加载基因映射字典...")
    gene_dict = load_gene_mapping_dict()
    ensembl_to_symbol = create_ensembl_to_symbol_mapping(gene_dict)
    
    if ensembl_to_symbol is None:
        logging.error("无法创建基因映射，退出分析")
        return
    
    # 加载数据
    logging.info("加载相似度数据和基因表达数据...")
    try:
        # 加载相似度数据
        similarity_df = pd.read_csv(similarity_file)
        
        # 加载基因表达数据
        adata = sc.read_h5ad(expression_file)
        
        logging.info("数据加载成功。")
    except Exception as e:
        logging.error(f"数据加载失败: {e}")
        return
    
    # 准备 AnnData 对象
    logging.info("准备 AnnData 对象用于差异表达分析...")
    
    try:
        # 筛选目标基因的数据
        target_gene_data = similarity_df[similarity_df['gene_name'] == target_gene].copy()
        
        if target_gene_data.empty:
            logging.error(f"未找到基因 '{target_gene}' 的数据")
            return
        
        # 根据相似度变化分组
        target_gene_data['response_group'] = target_gene_data['similarity_shift'].apply(
            lambda x: 'responder' if x > 0 else 'non_responder'
        )
        
        # 获取细胞ID
        cell_ids = target_gene_data['cell_id'].values
        
        # 筛选 AnnData 对象（使用细胞ID作为索引）
        adata_subset = adata[adata.obs.index.isin(cell_ids), :].copy()
        
        # 添加分组信息
        response_groups = target_gene_data.set_index('cell_id')['response_group']
        adata_subset.obs['response_group'] = response_groups.loc[adata_subset.obs.index].values
        
        logging.info(f"已筛选出 {adata_subset.n_obs} 个细胞用于分析。")
        logging.info(f"响应组: {sum(adata_subset.obs['response_group'] == 'responder')} 个细胞")
        logging.info(f"非响应组: {sum(adata_subset.obs['response_group'] == 'non_responder')} 个细胞")
        
    except Exception as e:
        logging.error(f"数据准备失败: {e}")
        return
    
    # 数据预处理
    logging.info("对数据进行标准化和对数化处理...")
    try:
        # 标准化数据
        sc.pp.normalize_total(adata_subset, target_sum=1e4)
        sc.pp.log1p(adata_subset)
        
        logging.info("数据处理完成。")
    except Exception as e:
        logging.error(f"数据处理失败: {e}")
        return
    
    # 差异表达分析
    logging.info("开始进行差异表达基因分析 (DEG)...")
    try:
        # 使用 scanpy 进行差异表达分析
        sc.tl.rank_genes_groups(
            adata_subset, 
            'response_group', 
            method='t-test',
            reference='non_responder'
        )
        
        logging.info("DEG 分析完成。")
    except Exception as e:
        logging.error(f"DEG 分析失败: {e}")
        return
    
    # 准备 GSEA 输入
    logging.info("准备 GSEA 的输入数据 (基因排序列表)...")
    try:
        # 获取差异表达结果
        result = adata_subset.uns['rank_genes_groups']
        
        # 提取 responder 组相对于 non_responder 组的结果
        genes = result['names']['responder']
        scores = result['scores']['responder']
        pvals = result['pvals']['responder']
        logfoldchanges = result['logfoldchanges']['responder']
        
        # 创建结果DataFrame
        deg_df = pd.DataFrame({
            'ensembl_id': genes,
            'score': scores,
            'pval': pvals,
            'logfoldchange': logfoldchanges
        })
        
        # 转换Ensembl ID为基因符号
        deg_df['gene_symbol'] = deg_df['ensembl_id'].map(ensembl_to_symbol)
        
        # 移除无法映射的基因
        deg_df = deg_df.dropna(subset=['gene_symbol'])
        
        # 对于有重复基因符号的情况，保留绝对score值最大的
        deg_df['abs_score'] = deg_df['score'].abs()
        deg_df = deg_df.sort_values('abs_score', ascending=False)
        deg_df = deg_df.drop_duplicates(subset=['gene_symbol'], keep='first')
        
        # 创建基因排序列表（用于GSEA）
        gene_rank_list = deg_df.set_index('gene_symbol')['score'].sort_values(ascending=False)
        
        logging.info(f"GSEA 输入列表准备完成，包含 {len(gene_rank_list)} 个基因。")
        logging.info(f"成功映射了 {len(deg_df)} / {len(genes)} 个基因")
        
    except Exception as e:
        logging.error(f"GSEA 输入准备失败: {e}")
        return
    
    # 执行 GSEA 分析
    logging.info("开始执行 GSEA 分析，使用本地基因集文件...")
    try:
        # 使用本地下载的 Hallmark 基因集
        gmt_file = 'h.all.v2023.2.Hs.symbols.gmt'
        
        gsea_result = gseapy.prerank(
            rnk=gene_rank_list,
            gene_sets=gmt_file,
            min_size=15,
            max_size=500,
            permutation_num=1000,
            outdir='results/gsea_results',
            seed=42,
            verbose=True
        )
        
        logging.info("GSEA 分析完成。")
        
        # 保存结果
        gsea_result.res2d.to_csv('results/gsea_hallmark_results.csv', index=False)
        logging.info("GSEA 结果已保存到 'results/gsea_hallmark_results.csv'")
        
        # 显示前10个最显著的通路
        significant_pathways = gsea_result.res2d[gsea_result.res2d['FDR q-val'] < 0.25]
        if not significant_pathways.empty:
            logging.info(f"发现 {len(significant_pathways)} 个显著富集的通路 (FDR < 0.25):")
            for _, pathway in significant_pathways.head(10).iterrows():
                logging.info(f"  {pathway['Term']}: NES={pathway['NES']:.3f}, FDR={pathway['FDR q-val']:.3f}")
        else:
            logging.info("未发现显著富集的通路 (FDR < 0.25)")
        
    except Exception as e:
        logging.error(f"分析过程中发生错误: {e}")
        return
    
    logging.info("通路富集分析脚本运行结束。")

if __name__ == "__main__":
    main()
