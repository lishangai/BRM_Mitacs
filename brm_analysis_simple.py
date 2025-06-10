
"""
BRM基因RNA-seq数据分析脚本（简化版）
基于illuminahiseq_rnaseqv2-RSEM_genes_normalized数据
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
import warnings
warnings.filterwarnings('ignore')

def load_and_preprocess_data(filepath):
    """
    加载和预处理RSEM标准化数据
    """
    print("正在加载数据...")
    try:
        # 读取Excel文件
        expr = pd.read_excel(filepath, index_col=0)
    except Exception as e:
        print(f"读取Excel文件失败: {e}")
        print("尝试使用pandas其他方法...")
        try:
            expr = pd.read_excel(filepath, engine='openpyxl', index_col=0)
        except:
            print("无法读取文件，请检查文件格式")
            return None, None
    
    print(f"原始数据维度: {expr.shape}")
    print(f"基因数: {expr.shape[0]}, 样本数: {expr.shape[1]}")
    
    # 检查BRM基因是否存在
    brm_genes = ['BRM', 'SMARCA2']
    brm_found = None
    for gene in brm_genes:
        if gene in expr.index:
            brm_found = gene
            print(f"找到BRM基因: {gene}")
            break
    
    if brm_found is None:
        print("警告: 未找到BRM或SMARCA2基因")
        # 寻找包含这些关键词的基因
        possible_brm = [g for g in expr.index if 'BRM' in str(g).upper() or 'SMARCA2' in str(g).upper()]
        if possible_brm:
            print(f"可能的BRM相关基因: {possible_brm[:5]}")
            brm_found = possible_brm[0]
        else:
            print("未找到任何BRM相关基因，将使用第一个基因作为示例")
            brm_found = expr.index[0]
    
    # 过滤低表达基因 (保留在≥10%样本中表达>1的基因)
    print("过滤低表达基因...")
    before_filter = expr.shape[0]
    expr_filtered = expr[(expr > 1).sum(axis=1) >= 0.1 * expr.shape[1]]
    after_filter = expr_filtered.shape[0]
    print(f"过滤后保留基因: {after_filter}/{before_filter} ({after_filter/before_filter*100:.1f}%)")
    
    return expr_filtered, brm_found

def create_brm_groups(expr, brm_gene):
    """
    根据BRM表达水平创建高低表达分组
    """
    print(f"基于{brm_gene}表达创建分组...")
    
    # log2转换
    brm_expr = np.log2(expr.loc[brm_gene] + 1)
    brm_median = np.median(brm_expr)
    
    # 创建分组
    groups = np.where(brm_expr > brm_median, "High", "Low")
    groups_series = pd.Series(groups, index=expr.columns)
    
    print(f"高表达组: {sum(groups == 'High')}个样本")
    print(f"低表达组: {sum(groups == 'Low')}个样本")
    print(f"BRM表达中位数: {brm_median:.3f}")
    
    return groups_series

def differential_expression_analysis(expr, groups):
    """
    差异表达分析
    """
    print("进行差异表达分析...")
    
    # 添加FDR校正函数
    def fdr_correction(pvals):
        """简单的FDR校正"""
        pvals = np.array(pvals)
        sorted_indices = np.argsort(pvals)
        sorted_pvals = pvals[sorted_indices]
        n = len(pvals)
        
        # Benjamini-Hochberg方法
        adjusted_pvals = np.empty(n)
        for i in range(n-1, -1, -1):
            if i == n-1:
                adjusted_pvals[sorted_indices[i]] = sorted_pvals[i]
            else:
                adjusted_pvals[sorted_indices[i]] = min(
                    sorted_pvals[i] * n / (i + 1),
                    adjusted_pvals[sorted_indices[i+1]]
                )
        
        return np.clip(adjusted_pvals, 0, 1)
    
    deg_results = []
    total_genes = len(expr.index)
    
    for i, gene in enumerate(expr.index):
        if i % 1000 == 0:
            print(f"处理进度: {i}/{total_genes}")
            
        # 获取高低表达组的数据
        high_samples = groups[groups == "High"].index
        low_samples = groups[groups == "Low"].index
        
        high = np.log2(expr.loc[gene, high_samples] + 1)
        low = np.log2(expr.loc[gene, low_samples] + 1)
        
        # Mann-Whitney U检验
        try:
            statistic, pval = mannwhitneyu(high, low, alternative="two-sided")
            logfc = np.mean(high) - np.mean(low)
            deg_results.append([gene, logfc, pval, np.mean(high), np.mean(low)])
        except:
            deg_results.append([gene, 0, 1, 0, 0])
    
    # 创建结果DataFrame
    deg_df = pd.DataFrame(deg_results, columns=["gene", "log2FC", "pval", "mean_high", "mean_low"])
    
    # FDR校正
    deg_df["padj"] = fdr_correction(deg_df["pval"])
    
    # 筛选显著差异基因
    significant_genes = deg_df[deg_df["padj"] < 0.05].copy()
    significant_genes = significant_genes.reindex(significant_genes["log2FC"].abs().sort_values(ascending=False).index)
    
    print(f"显著差异基因数: {len(significant_genes)}")
    print(f"上调基因数: {sum(significant_genes['log2FC'] > 0)}")
    print(f"下调基因数: {sum(significant_genes['log2FC'] < 0)}")
    
    return deg_df, significant_genes

def correlation_analysis(expr, brm_gene):
    """
    与BRM基因的相关性分析
    """
    print(f"计算与{brm_gene}的相关性...")
    
    # 添加FDR校正函数
    def fdr_correction(pvals):
        """简单的FDR校正"""
        pvals = np.array(pvals)
        sorted_indices = np.argsort(pvals)
        sorted_pvals = pvals[sorted_indices]
        n = len(pvals)
        
        # Benjamini-Hochberg方法
        adjusted_pvals = np.empty(n)
        for i in range(n-1, -1, -1):
            if i == n-1:
                adjusted_pvals[sorted_indices[i]] = sorted_pvals[i]
            else:
                adjusted_pvals[sorted_indices[i]] = min(
                    sorted_pvals[i] * n / (i + 1),
                    adjusted_pvals[sorted_indices[i+1]]
                )
        
        return np.clip(adjusted_pvals, 0, 1)
    
    corr_results = []
    brm_expr = expr.loc[brm_gene]
    
    for gene in expr.index:
        try:
            rho, pval = spearmanr(brm_expr, expr.loc[gene])
            corr_results.append([gene, rho, pval])
        except:
            corr_results.append([gene, 0, 1])
    
    corr_df = pd.DataFrame(corr_results, columns=["gene", "rho", "pval"])
    corr_df["padj"] = fdr_correction(corr_df["pval"])
    
    # 筛选显著相关基因
    significant_corr = corr_df[(corr_df["padj"] < 0.05) & (abs(corr_df["rho"]) > 0.3)]
    
    print(f"显著正相关基因数: {sum((significant_corr['rho'] > 0.3) & (significant_corr['padj'] < 0.05))}")
    print(f"显著负相关基因数: {sum((significant_corr['rho'] < -0.3) & (significant_corr['padj'] < 0.05))}")
    
    return corr_df, significant_corr

def save_results(deg_df, significant_genes, corr_df, significant_corr):
    """
    保存分析结果
    """
    print("保存分析结果...")
    
    # 保存差异表达结果
    deg_df.to_csv('differential_expression_results.csv', index=False)
    significant_genes.to_csv('significant_DEGs.csv', index=False)
    
    # 保存相关性分析结果
    corr_df.to_csv('correlation_analysis_results.csv', index=False)
    significant_corr.to_csv('significant_correlations.csv', index=False)
    
    print("结果已保存:")
    print("- differential_expression_results.csv: 所有基因差异表达结果")
    print("- significant_DEGs.csv: 显著差异表达基因")
    print("- correlation_analysis_results.csv: 所有基因相关性结果")
    print("- significant_correlations.csv: 显著相关基因")

def print_summary(brm_gene, expr, significant_genes, significant_corr, deg_df):
    """
    打印分析总结
    """
    print("\n=== 分析总结 ===")
    print(f"分析基因: {brm_gene}")
    print(f"总基因数: {expr.shape[0]}")
    print(f"总样本数: {expr.shape[1]}")
    print(f"显著差异基因数: {len(significant_genes)}")
    print(f"显著相关基因数: {len(significant_corr)}")
    
    if len(significant_genes) > 0:
        print(f"\nTop 10 上调基因:")
        top_up = significant_genes[significant_genes['log2FC'] > 0].head(10)
        for _, row in top_up.iterrows():
            print(f"  {row['gene']}: log2FC={row['log2FC']:.3f}, padj={row['padj']:.2e}")
        
        print(f"\nTop 10 下调基因:")
        top_down = significant_genes[significant_genes['log2FC'] < 0].head(10)
        for _, row in top_down.iterrows():
            print(f"  {row['gene']}: log2FC={row['log2FC']:.3f}, padj={row['padj']:.2e}")
    
    # 统计信息
    print(f"\n=== 统计信息 ===")
    print(f"总p值 < 0.05的基因数: {sum(deg_df['pval'] < 0.05)}")
    print(f"FDR校正后 < 0.05的基因数: {sum(deg_df['padj'] < 0.05)}")
    print(f"log2FC > 1的基因数: {sum(abs(deg_df['log2FC']) > 1)}")
    print(f"log2FC > 2的基因数: {sum(abs(deg_df['log2FC']) > 2)}")

def main():
    """
    主分析流程
    """
    print("=== BRM基因RNA-seq数据分析 ===\n")
    
    # 1. 数据加载和预处理
    expr, brm_gene = load_and_preprocess_data("illuminahiseq_rnaseqv2-RSEM_genes_normalized (MD5).xlsx")
    
    if expr is None:
        print("数据加载失败，退出分析")
        return
    
    # 2. 创建BRM高低表达分组
    groups = create_brm_groups(expr, brm_gene)
    
    # 3. 差异表达分析
    deg_df, significant_genes = differential_expression_analysis(expr, groups)
    
    # 4. 相关性分析
    corr_df, significant_corr = correlation_analysis(expr, brm_gene)
    
    # 5. 保存结果
    save_results(deg_df, significant_genes, corr_df, significant_corr)
    
    # 6. 输出分析总结
    print_summary(brm_gene, expr, significant_genes, significant_corr, deg_df)

if __name__ == "__main__":
    main() 