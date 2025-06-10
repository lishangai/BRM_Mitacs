
"""
BRM基因RNA-seq数据分析脚本（支持CSV格式）
基于illuminahiseq_rnaseqv2-RSEM_genes_normalized数据
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
import warnings
import os
warnings.filterwarnings('ignore')

def load_and_preprocess_data(filepath):
    """
    加载和预处理RSEM标准化数据
    """
    print("正在加载数据...")
    
    # 检查文件扩展名
    if filepath.endswith('.csv'):
        try:
            expr = pd.read_csv(filepath, index_col=0)
            print("成功读取CSV文件")
        except Exception as e:
            print(f"读取CSV文件失败: {e}")
            return None, None
    elif filepath.endswith('.xlsx') or filepath.endswith('.xls'):
        try:
            expr = pd.read_excel(filepath, index_col=0)
            print("成功读取Excel文件")
        except Exception as e:
            print(f"读取Excel文件失败: {e}")
            # 尝试其他方法
            try:
                expr = pd.read_excel(filepath, engine='openpyxl', index_col=0)
                print("使用openpyxl引擎成功读取Excel文件")
            except:
                print("Excel文件读取失败，请将文件转换为CSV格式")
                return None, None
    else:
        print("不支持的文件格式，请使用CSV或Excel文件")
        return None, None
    
    print(f"原始数据维度: {expr.shape}")
    print(f"基因数: {expr.shape[0]}, 样本数: {expr.shape[1]}")
    
    # 显示数据的前几行和前几列
    print(f"\n数据预览:")
    print(expr.iloc[:5, :5])
    
    # 数据清理：确保所有列都是数值型
    print("清理数据，转换为数值型...")
    
    # 删除第一行如果它包含非数值信息
    if 'Hybridization REF' in expr.index:
        expr = expr.drop('Hybridization REF')
        print("删除了Hybridization REF行")
    
    # 删除可能的描述行
    non_numeric_rows = []
    for idx in expr.index:
        if idx in ['gene_id', 'normalized_count', 'Hybridization REF']:
            non_numeric_rows.append(idx)
    
    if non_numeric_rows:
        expr = expr.drop(non_numeric_rows)
        print(f"删除了非数值行: {non_numeric_rows}")
    
    # 将所有数据转换为数值型
    try:
        expr = expr.apply(pd.to_numeric, errors='coerce')
        # 删除包含NaN的行
        before_na = expr.shape[0]
        expr = expr.dropna()
        after_na = expr.shape[0]
        if before_na != after_na:
            print(f"删除了包含缺失值的行: {before_na - after_na}行")
    except Exception as e:
        print(f"数据转换警告: {e}")
    
    print(f"清理后数据维度: {expr.shape}")
    
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
            # 优先选择SMARCA2
            smarca2_genes = [g for g in possible_brm if 'SMARCA2' in str(g).upper()]
            if smarca2_genes:
                brm_found = smarca2_genes[0]
                print(f"选择SMARCA2基因: {brm_found}")
            else:
                brm_found = possible_brm[0]
                print(f"未找到SMARCA2，选择: {brm_found}")
        else:
            print("未找到任何BRM相关基因")
            # 显示一些基因名例子
            print(f"基因名示例: {list(expr.index[:10])}")
            print("请确认BRM基因的正确名称，或手动选择一个基因进行分析")
            brm_found = expr.index[0]  # 使用第一个基因作为示例
            print(f"使用 {brm_found} 作为示例基因进行分析")
    
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
    print(f"BRM表达范围: {brm_expr.min():.3f} - {brm_expr.max():.3f}")
    
    return groups_series

def differential_expression_analysis(expr, groups):
    """
    差异表达分析
    """
    print("进行差异表达分析...")
    
    # FDR校正函数
    def fdr_correction(pvals):
        """Benjamini-Hochberg FDR校正"""
        pvals = np.array(pvals)
        sorted_indices = np.argsort(pvals)
        sorted_pvals = pvals[sorted_indices]
        n = len(pvals)
        
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
    
    # FDR校正函数
    def fdr_correction(pvals):
        """Benjamini-Hochberg FDR校正"""
        pvals = np.array(pvals)
        sorted_indices = np.argsort(pvals)
        sorted_pvals = pvals[sorted_indices]
        n = len(pvals)
        
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
    print(f"|log2FC| > 1的基因数: {sum(abs(deg_df['log2FC']) > 1)}")
    print(f"|log2FC| > 2的基因数: {sum(abs(deg_df['log2FC']) > 2)}")

def main():
    """
    主分析流程
    """
    print("=== BRM基因RNA-seq数据分析 ===\n")
    
    # 检查可能的数据文件
    possible_files = [
        "illuminahiseq_rnaseqv2-RSEM_genes_normalized.csv",
        "illuminahiseq_rnaseqv2-RSEM_genes_normalized (MD5).csv",
        "illuminahiseq_rnaseqv2-RSEM_genes_normalized (MD5).xlsx",
        "illuminahiseq_rnaseqv2-RSEM_genes_normalized.xlsx"
    ]
    
    filepath = None
    for file in possible_files:
        if os.path.exists(file):
            filepath = file
            print(f"找到数据文件: {filepath}")
            break
    
    if filepath is None:
        print("未找到数据文件，请确保以下文件之一存在:")
        for file in possible_files:
            print(f"  - {file}")
        return
    
    # 1. 数据加载和预处理
    expr, brm_gene = load_and_preprocess_data(filepath)
    
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