
"""
BRM基因RNA-seq数据分析脚本
基于illuminahiseq_rnaseqv2-RSEM_genes_normalized数据
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import fdrcorrection
import warnings
warnings.filterwarnings('ignore')

# 设置图形样式
plt.style.use('default')
sns.set_palette("husl")

def load_and_preprocess_data(filepath):
    """
    加载和预处理RSEM标准化数据
    """
    print("正在加载数据...")
    # 读取Excel文件
    expr = pd.read_excel(filepath, index_col=0)
    
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
    
    # 可视化BRM表达分布
    plt.figure(figsize=(10, 6))
    plt.subplot(1, 2, 1)
    plt.hist(brm_expr, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(brm_median, color='red', linestyle='--', label=f'中位数: {brm_median:.2f}')
    plt.xlabel(f'{brm_gene} Expression (log2(TPM+1))')
    plt.ylabel('样本数')
    plt.title(f'{brm_gene} 表达分布')
    plt.legend()
    
    plt.subplot(1, 2, 2)
    box_data = [brm_expr[groups == 'Low'], brm_expr[groups == 'High']]
    plt.boxplot(box_data, labels=['Low', 'High'])
    plt.ylabel(f'{brm_gene} Expression (log2(TPM+1))')
    plt.title(f'{brm_gene} 分组表达水平')
    
    plt.tight_layout()
    plt.savefig('brm_expression_distribution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return groups_series

def differential_expression_analysis(expr, groups):
    """
    差异表达分析
    """
    print("进行差异表达分析...")
    
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
    deg_df["padj"] = fdrcorrection(deg_df["pval"])[1]
    
    # 筛选显著差异基因
    significant_genes = deg_df[deg_df["padj"] < 0.05].copy()
    significant_genes = significant_genes.sort_values("log2FC", key=abs, ascending=False)
    
    print(f"显著差异基因数: {len(significant_genes)}")
    print(f"上调基因数: {sum(significant_genes['log2FC'] > 0)}")
    print(f"下调基因数: {sum(significant_genes['log2FC'] < 0)}")
    
    return deg_df, significant_genes

def create_volcano_plot(deg_df):
    """
    创建火山图
    """
    print("绘制火山图...")
    
    plt.figure(figsize=(10, 8))
    
    # 计算-log10(padj)
    deg_df['-log10_padj'] = -np.log10(deg_df['padj'] + 1e-300)  # 避免log(0)
    
    # 设置颜色
    colors = []
    for _, row in deg_df.iterrows():
        if row['padj'] < 0.05 and abs(row['log2FC']) > 1:
            if row['log2FC'] > 0:
                colors.append('red')  # 上调
            else:
                colors.append('blue')  # 下调
        else:
            colors.append('gray')  # 不显著
    
    # 绘制散点图
    plt.scatter(deg_df['log2FC'], deg_df['-log10_padj'], 
               c=colors, alpha=0.6, s=10)
    
    # 添加阈值线
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=1, color='black', linestyle='--', alpha=0.5)
    plt.axvline(x=-1, color='black', linestyle='--', alpha=0.5)
    
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(adjusted p-value)')
    plt.title('Volcano Plot: BRM-High vs BRM-Low')
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='red', label='上调显著'),
                      Patch(facecolor='blue', label='下调显著'),
                      Patch(facecolor='gray', label='不显著')]
    plt.legend(handles=legend_elements)
    
    plt.tight_layout()
    plt.savefig('volcano_plot.png', dpi=300, bbox_inches='tight')
    plt.show()

def create_deg_heatmap(expr, significant_genes, groups, top_n=50):
    """
    创建Top DEGs热图
    """
    print(f"绘制Top {top_n} 差异基因热图...")
    
    if len(significant_genes) == 0:
        print("没有显著差异基因，跳过热图绘制")
        return
    
    # 选择top基因
    top_genes = significant_genes.head(top_n)["gene"].tolist()
    
    # 准备数据
    heatmap_data = np.log2(expr.loc[top_genes] + 1)
    
    # 准备列颜色
    col_colors = groups.map({"High": "red", "Low": "blue"})
    
    # 绘制聚类热图
    try:
        g = sns.clustermap(heatmap_data, 
                          col_colors=col_colors,
                          figsize=(12, 10),
                          cmap='RdBu_r',
                          center=0,
                          xticklabels=False,
                          yticklabels=True)
        
        g.savefig('deg_heatmap.png', dpi=300, bbox_inches='tight')
        plt.show()
    except Exception as e:
        print(f"热图绘制失败: {e}")

def correlation_analysis(expr, brm_gene):
    """
    与BRM基因的相关性分析
    """
    print(f"计算与{brm_gene}的相关性...")
    
    corr_results = []
    brm_expr = expr.loc[brm_gene]
    
    for gene in expr.index:
        try:
            rho, pval = spearmanr(brm_expr, expr.loc[gene])
            corr_results.append([gene, rho, pval])
        except:
            corr_results.append([gene, 0, 1])
    
    corr_df = pd.DataFrame(corr_results, columns=["gene", "rho", "pval"])
    corr_df["padj"] = fdrcorrection(corr_df["pval"])[1]
    
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

def main():
    """
    主分析流程
    """
    print("=== BRM基因RNA-seq数据分析 ===\n")
    
    # 1. 数据加载和预处理
    expr, brm_gene = load_and_preprocess_data("illuminahiseq_rnaseqv2-RSEM_genes_normalized (MD5).xlsx")
    
    # 2. 创建BRM高低表达分组
    groups = create_brm_groups(expr, brm_gene)
    
    # 3. 差异表达分析
    deg_df, significant_genes = differential_expression_analysis(expr, groups)
    
    # 4. 可视化结果
    create_volcano_plot(deg_df)
    create_deg_heatmap(expr, significant_genes, groups)
    
    # 5. 相关性分析
    corr_df, significant_corr = correlation_analysis(expr, brm_gene)
    
    # 6. 保存结果
    save_results(deg_df, significant_genes, corr_df, significant_corr)
    
    # 7. 输出分析总结
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
            print(f"  {row['gene']}: log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")
        
        print(f"\nTop 10 下调基因:")
        top_down = significant_genes[significant_genes['log2FC'] < 0].head(10)
        for _, row in top_down.iterrows():
            print(f"  {row['gene']}: log2FC={row['log2FC']:.2f}, padj={row['padj']:.2e}")

if __name__ == "__main__":
    main() 