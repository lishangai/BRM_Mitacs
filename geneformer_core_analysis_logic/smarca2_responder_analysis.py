import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import scanpy as sc
import pandas as pd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

def setup_logging_and_dirs(gene_name):
    """设置日志和输出目录"""
    script_name = Path(__file__).stem
    project_root = Path(__file__).resolve().parent.parent
    base_results_dir = project_root / 'results' / script_name
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    gene_specific_dir = base_results_dir / f"{gene_name}_{timestamp}"
    gene_specific_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = gene_specific_dir / "analysis_log.txt"
    
    def log_func(message):
        timestamp_log = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
        log_message = f"{timestamp_log} {message}"
        print(log_message)
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(log_message + '\n')
            
    return log_func, gene_specific_dir

def load_data(log_func):
    """加载所需数据"""
    log_func("Loading data...")
    project_root = Path(__file__).resolve().parent.parent
    
    # 加载相似度变化结果
    shift_results_path = project_root / 'results' / 'refractory_to_sensitive_shift_analysis' / 'similarity_shift_results.csv'
    if not shift_results_path.exists():
        raise FileNotFoundError(f"Similarity shift data not found at: {shift_results_path}")
    shift_df = pd.read_csv(shift_results_path)
    log_func(f"Loaded similarity shift data. Shape: {shift_df.shape}")

    # 加载原始 AnnData 对象
    adata_path = project_root / "data" / "processed" / "geneformer_input.h5ad"
    if not adata_path.exists():
        raise FileNotFoundError(f"AnnData object not found at: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    log_func(f"Loaded AnnData object. Shape: {adata.shape}")
    
    return shift_df, adata

def identify_and_group_cells(shift_df, adata, gene_name, log_func):
    """识别并分组特定基因的细胞"""
    log_func(f"Identifying and grouping cells for gene: {gene_name}")
    
    # 筛选特定基因的细胞
    gene_df = shift_df[shift_df['gene_name'] == gene_name].copy()
    if gene_df.empty:
        raise ValueError(f"No data found for gene '{gene_name}' in the similarity shift results.")
        
    log_func(f"Found {len(gene_df)} cells analyzed for {gene_name} knockout.")
    
    # 分组
    gene_df['group'] = np.where(gene_df['similarity_shift'] > 0, 'Responder', 'Non-responder')
    
    # 统计分组情况
    responder_count = (gene_df['group'] == 'Responder').sum()
    non_responder_count = (gene_df['group'] == 'Non-responder').sum()
    log_func(f"Responders (Similarity Shift > 0): {responder_count} cells")
    log_func(f"Non-responders (Similarity Shift <= 0): {non_responder_count} cells")

    # 在 AnnData 对象中添加分组信息
    cell_barcodes = gene_df['cell_id'].tolist()
    adata_subset = adata[adata.obs.index.isin(cell_barcodes)].copy()
    
    group_map = gene_df.set_index('cell_id')['group']
    adata_subset.obs['responder_status'] = adata_subset.obs.index.map(group_map)
    
    return adata_subset

def perform_deg_analysis(adata_subset, log_func, results_dir):
    """执行差异表达基因分析"""
    log_func("Performing Differential Gene Expression (DEG) analysis...")
    
    # 确保有两组可供比较
    if len(adata_subset.obs['responder_status'].unique()) < 2:
        log_func("Warning: Only one group present. Skipping DEG analysis.")
        return None, None
        
    # 预处理
    sc.pp.normalize_total(adata_subset, target_sum=1e4)
    sc.pp.log1p(adata_subset)
    
    # DEG 分析
    log_func("Running rank_genes_groups (t-test)...")
    sc.tl.rank_genes_groups(adata_subset, 'responder_status', method='t-test_overestim_var')
    
    # 跳过可视化以避免 NumPy 2.0 兼容性问题
    log_func("Skipping built-in visualization due to NumPy 2.0 compatibility issues.")
    
    # 提取结果
    result = adata_subset.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    deg_results = {}
    for group in groups:
        deg_results[group] = pd.DataFrame({
            'names': result['names'][group],
            'scores': result['scores'][group],
            'pvals': result['pvals'][group],
            'pvals_adj': result['pvals_adj'][group],
            'logfoldchanges': result['logfoldchanges'][group],
        })
        
    responder_degs = deg_results.get('Responder')
    if responder_degs is None:
        log_func("Could not find 'Responder' group in DEG results.")
        return None, None

    # 保存完整的 DEG 列表
    deg_list_path = results_dir / 'responder_vs_non-responder_deg.csv'
    responder_degs.to_csv(deg_list_path, index=False)
    log_func(f"Full DEG list saved to: {deg_list_path}")
    
    return responder_degs, adata_subset

def create_volcano_plot(deg_df, log_func, results_dir):
    """创建火山图"""
    if deg_df is None:
        return
        
    log_func("Creating volcano plot...")
    
    deg_df['minus_log10_p_adj'] = -np.log10(deg_df['pvals_adj'])
    
    # 增加一个阈值，避免 -log10(0) 变成无穷大
    deg_df.replace([np.inf, -np.inf], 350, inplace=True) # 350是一个常用的上限

    plt.figure(figsize=(10, 8))
    
    # 默认颜色
    plt.scatter(deg_df['logfoldchanges'], deg_df['minus_log10_p_adj'], 
                alpha=0.4, s=10, c='grey')

    # 标记显著上调和下调的基因
    up_regulated = deg_df[(deg_df['logfoldchanges'] > 1) & (deg_df['pvals_adj'] < 0.05)]
    down_regulated = deg_df[(deg_df['logfoldchanges'] < -1) & (deg_df['pvals_adj'] < 0.05)]

    plt.scatter(up_regulated['logfoldchanges'], up_regulated['minus_log10_p_adj'], 
                alpha=0.6, s=20, c='red', label='Up-regulated in Responders')
    plt.scatter(down_regulated['logfoldchanges'], down_regulated['minus_log10_p_adj'], 
                alpha=0.6, s=20, c='blue', label='Down-regulated in Responders')

    # 添加基因标签
    genes_to_label = pd.concat([
        up_regulated.nlargest(10, 'minus_log10_p_adj'),
        down_regulated.nlargest(10, 'minus_log10_p_adj')
    ])
    
    for i, row in genes_to_label.iterrows():
        plt.text(row['logfoldchanges'], row['minus_log10_p_adj'], row['names'],
                 fontsize=9, ha='center', va='bottom')

    plt.title('Responders vs. Non-responders', fontsize=16, fontweight='bold')
    plt.xlabel('Log2 Fold Change', fontsize=14)
    plt.ylabel('-log10(Adjusted p-value)', fontsize=14)
    plt.axvline(x=1, color='grey', linestyle='--', lw=1)
    plt.axvline(x=-1, color='grey', linestyle='--', lw=1)
    plt.axhline(y=-np.log10(0.05), color='grey', linestyle='--', lw=1)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.2)
    
    volcano_path = results_dir / 'volcano_plot.png'
    plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
    plt.close()
    log_func(f"Volcano plot saved to: {volcano_path}")

def generate_summary(deg_df, log_func, gene_name):
    """生成并打印摘要"""
    if deg_df is None:
        log_func("No DEG results to summarize.")
        return
        
    log_func("\n" + "="*30)
    log_func(f"Summary for {gene_name} Responder Analysis")
    log_func("="*30)
    
    up_regulated = deg_df[(deg_df['logfoldchanges'] > 1) & (deg_df['pvals_adj'] < 0.05)]
    down_regulated = deg_df[(deg_df['logfoldchanges'] < -1) & (deg_df['pvals_adj'] < 0.05)]
    
    log_func(f"Found {len(up_regulated)} significantly UP-regulated genes in Responders.")
    log_func(f"Found {len(down_regulated)} significantly DOWN-regulated genes in Responders.")
    
    log_func("\nTop 10 UP-regulated genes in Responders (potential synergistic targets for knockout):")
    log_func("-" * 80)
    log_func(up_regulated.head(10).to_string())
    
    log_func("\nTop 10 DOWN-regulated genes in Responders (genes whose absence may enable positive response):")
    log_func("-" * 80)
    log_func(down_regulated.head(10).to_string())
    log_func("\n")

def main():
    GENE_OF_INTEREST = "SMARCA2"
    log_func, results_dir = setup_logging_and_dirs(GENE_OF_INTEREST)
    
    try:
        shift_df, adata = load_data(log_func)
        adata_subset = identify_and_group_cells(shift_df, adata, GENE_OF_INTEREST, log_func)
        deg_results, adata_with_degs = perform_deg_analysis(adata_subset, log_func, results_dir)
        create_volcano_plot(deg_results, log_func, results_dir)
        generate_summary(deg_results, log_func, GENE_OF_INTEREST)
        
        log_func("Analysis complete.")
    except (FileNotFoundError, ValueError) as e:
        log_func(f"Error: {e}")
    except Exception as e:
        log_func(f"An unexpected error occurred: {e}")
        import traceback
        log_func(traceback.format_exc())

if __name__ == "__main__":
    main()
