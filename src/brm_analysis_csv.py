 
"""
BRM基因RNA-seq数据分析脚本（支持CSV格式）
基于illuminahiseq_rnaseqv2-RSEM_genes_normalized数据
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, spearmanr
import warnings
import os
from pathlib import Path
import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection

warnings.filterwarnings('ignore')

# Define Project Root for robust path management
PROJECT_ROOT = Path(__file__).resolve().parent.parent

class Tee:
    """A helper class to redirect print output to both console and a file."""
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()  # Ensure output is written immediately
    def flush(self):
        for f in self.files:
            f.flush()

def plot_brm_distribution(brm_expr_values, brm_gene, output_dir):
    """
    绘制BRM基因表达分布的直方图，并标注统计信息
    """
    print(f"为 {brm_gene} 表达分布创建可视化图表...")
    
    # 计算统计数据
    mean_val = np.mean(brm_expr_values)
    median_val = np.median(brm_expr_values)
    q1 = np.percentile(brm_expr_values, 25)
    q3 = np.percentile(brm_expr_values, 75)
    
    plt.figure(figsize=(12, 8))
    
    # 绘制直方图和KDE曲线
    sns.histplot(brm_expr_values, kde=False, bins=30, color='skyblue', label='Sample Distribution', stat="density")
    
    # 拟合正态分布并绘制
    mu, std = norm.fit(brm_expr_values)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2, label='Normal dist. fit')
    
    # 标注统计线
    plt.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}')
    plt.axvline(median_val, color='green', linestyle='-', linewidth=2, label=f'Median: {median_val:.2f}')
    plt.axvline(q1, color='orange', linestyle=':', linewidth=2, label=f'1st Quartile (Q1): {q1:.2f}')
    plt.axvline(q3, color='purple', linestyle=':', linewidth=2, label=f'3rd Quartile (Q3): {q3:.2f}')
    
    # Sanitize the gene name to be safe for filenames and labels
    if '|' in brm_gene:
        safe_brm_gene_name = brm_gene.split('|')[0]
    else:
        safe_brm_gene_name = brm_gene.replace('|', '_').replace('?', '_')

    # 美化图表
    plt.title(f'Distribution of {safe_brm_gene_name} Expression', fontsize=16)
    plt.xlabel(f'Log2-Transformed Expression of {safe_brm_gene_name} (log2(value + 1))', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    # 保存图表
    plot_path = output_dir / f"{safe_brm_gene_name}_expression_distribution.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"分布图已保存到: {plot_path}")
    plt.close()

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
            except Exception as e2:
                print(f"Excel文件读取失败: {e2}")
                return None, None
    else:
        print("不支持的文件格式，请使用CSV或Excel文件")
        return None, None
    
    print(f"原始数据维度: {expr.shape}")
    print(f"基因数: {expr.shape[0]}, 样本数: {expr.shape[1]}")
    
    # 显示数据的前几行和前几列
    print(f"\n数据预览 (加载后):")
    print(expr.head())
    
    # 优先处理已知的非数值行
    # 删除第一行如果它包含非数值信息
    if 'Hybridization REF' in expr.index:
        print("\n检测到 'Hybridization REF' 行，准备删除...")
        print("删除前的索引预览:", expr.index[:5])
        expr = expr.drop('Hybridization REF')
        print("已删除 'Hybridization REF' 行。")
        print("删除后的索引预览:", expr.index[:5])
    
    # 删除可能的描述行
    non_numeric_rows = []
    for idx in expr.index:
        if idx in ['gene_id', 'normalized_count', 'Hybridization REF']:
            non_numeric_rows.append(idx)
    
    if non_numeric_rows:
        expr = expr.drop(non_numeric_rows)
        print(f"删除了非数值行: {non_numeric_rows}")
        
    print("\n打印纯数据预览 (无索引, 5x5):")
    print(expr.values[:5, :5])
        
    print("\n打印行索引预览 (前5个):")
    print(expr.index[:5])
    
    print("\n打印列索引预览 (前5个):")
    print(expr.columns[:5])
    
    # 数据清理：确保所有列都是数值型
    print("清理数据，转换为数值型...")
    
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

    
    # 新增：过滤0值占比过高的列和行
    print("过滤0值占比 > 50%的列（样本）...")
    before_filter_cols = expr.shape[1]
    col_zero_prop = (expr == 0).mean(axis=0)
    cols_to_keep = col_zero_prop[col_zero_prop <= 0.5].index
    expr = expr[cols_to_keep]
    after_filter_cols = expr.shape[1]
    print(f"删除了 {before_filter_cols - after_filter_cols} 个0值占比 > 50%的列")

    print("过滤0值占比 > 50%的行（基因）...")
    before_filter_rows = expr.shape[0]
    row_zero_prop = (expr == 0).mean(axis=1)
    rows_to_keep = row_zero_prop[row_zero_prop <= 0.5].index
    expr = expr.loc[rows_to_keep]
    after_filter_rows = expr.shape[0]
    print(f"删除了 {before_filter_rows - after_filter_rows} 个0值占比 > 50%的行")
        
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
    
    return expr, brm_found

def create_brm_groups(brm_expr, brm_gene):
    """
    根据BRM表达水平创建高低表达分组 (低表达组：最低40%，高表达组：最高40%)
    """
    print(f"基于 {brm_gene} 表达水平创建分组...")
    
    # 计算分位数：低40%和高40%
    q_low = np.percentile(brm_expr, 40)   # 40th percentile
    q_high = np.percentile(brm_expr, 60)  # 60th percentile
    
    # 创建分组
    def assign_group(value):
        if value <= q_low:
            return "Low"
        elif value >= q_high:
            return "High"
        else:
            return "Middle"

    groups = brm_expr.apply(assign_group)
    
    print(f"低表达组 (<= 40%, {q_low:.3f}): {sum(groups == 'Low')}个样本")
    print(f"高表达组 (>= 60%, {q_high:.3f}): {sum(groups == 'High')}个样本")
    print(f"中间组: {sum(groups == 'Middle')}个样本")
    
    return groups

def differential_expression_analysis(expr, groups, brm_gene):
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
    
    # FDR校正 (使用statsmodels的公认函数)
    _, deg_df["padj"] = fdrcorrection(deg_df["pval"], alpha=0.05)
    
    # 筛选显著差异基因
    significant_genes = deg_df[deg_df["padj"] < 0.05].copy()
    
    # 去除目标基因本身
    significant_genes = significant_genes[significant_genes['gene'] != brm_gene]
    
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
    
    corr_results = []
    brm_expr = expr.loc[brm_gene]
    
    for gene in expr.index:
        try:
            rho, pval = spearmanr(brm_expr, expr.loc[gene])
            corr_results.append([gene, rho, pval])
        except:
            corr_results.append([gene, 0, 1])
    
    corr_df = pd.DataFrame(corr_results, columns=["gene", "rho", "pval"])
    
    # FDR校正 (使用statsmodels的公认函数)
    _, corr_df["padj"] = fdrcorrection(corr_df["pval"], alpha=0.05)
    
    # 筛选显著相关基因
    significant_corr = corr_df[(corr_df["padj"] < 0.05) & (abs(corr_df["rho"]) > 0.3)].copy()
    
    # 去除目标基因本身
    significant_corr = significant_corr[significant_corr['gene'] != brm_gene]
    
    print(f"显著正相关基因数: {sum((significant_corr['rho'] > 0.3) & (significant_corr['padj'] < 0.05))}")
    print(f"显著负相关基因数: {sum((significant_corr['rho'] < -0.3) & (significant_corr['padj'] < 0.05))}")
    
    return corr_df, significant_corr

def save_results(deg_df, significant_genes, corr_df, significant_corr, output_dir):
    """
    保存分析结果
    """
    print(f"保存分析结果到: {output_dir}")
    
    # 确保输出目录存在
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 保存差异表达结果
    deg_df.to_csv(output_dir / 'differential_expression_results.csv', index=False)
    significant_genes.to_csv(output_dir / 'significant_DEGs.csv', index=False)
    
    # 保存相关性分析结果
    corr_df.to_csv(output_dir / 'correlation_analysis_results.csv', index=False)
    significant_corr.to_csv(output_dir / 'significant_correlations.csv', index=False)
    
    print("结果已保存:")
    print(f"- {output_dir / 'differential_expression_results.csv'}")
    print(f"- {output_dir / 'significant_DEGs.csv'}")
    print(f"- {output_dir / 'correlation_analysis_results.csv'}")
    print(f"- {output_dir / 'significant_correlations.csv'}")

def print_summary(brm_gene, expr, significant_genes, significant_corr, deg_df, corr_df):
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
        print(f"\nTop 15 显著上调基因 (padj < 0.05):")
        top_up = significant_genes[significant_genes['log2FC'] > 0].head(15)
        for _, row in top_up.iterrows():
            print(f"  {row['gene']}: log2FC={row['log2FC']:.3f}, padj={row['padj']:.2e}")
        
        print(f"\nTop 15 显著下调基因 (padj < 0.05):")
        top_down = significant_genes[significant_genes['log2FC'] < 0].head(15)
        for _, row in top_down.iterrows():
            print(f"  {row['gene']}: log2FC={row['log2FC']:.3f}, padj={row['padj']:.2e}")
    else:
        print(f"\n未找到显著差异表达基因 (padj < 0.05)")
    
    # 添加相关性分析结果 (仅显示显著性基因 padj < 0.05)
    if len(significant_corr) > 0:
        print(f"\nTop 15 显著正相关基因 (padj < 0.05, rho > 0.3):")
        # 按相关系数降序排列，取前10个正相关显著基因
        top_pos_corr = significant_corr[significant_corr['rho'] > 0].sort_values('rho', ascending=False).head(15)
        for _, row in top_pos_corr.iterrows():
            print(f"  {row['gene']}: rho={row['rho']:.3f}, padj={row['padj']:.2e}")
        
        print(f"\nTop 15 显著负相关基因 (padj < 0.05, rho < -0.3):")
        # 按相关系数升序排列，取前10个负相关显著基因
        top_neg_corr = significant_corr[significant_corr['rho'] < 0].sort_values('rho', ascending=True).head(15)
        for _, row in top_neg_corr.iterrows():
            print(f"  {row['gene']}: rho={row['rho']:.3f}, padj={row['padj']:.2e}")
    else:
        print(f"\n未找到显著相关基因 (padj < 0.05, |rho| > 0.3)")
    
    # 统计信息
    print(f"\n=== 统计信息 ===")
    print(f"总p值 < 0.05的基因数: {sum(deg_df['pval'] < 0.05)}")
    print(f"FDR校正后 < 0.05的基因数: {sum(deg_df['padj'] < 0.05)}")
    print(f"|log2FC| > 1的基因数: {sum(abs(deg_df['log2FC']) > 1)}")
    print(f"|log2FC| > 2的基因数: {sum(abs(deg_df['log2FC']) > 2)}")
    
    # 相关性统计信息
    print(f"总相关性p值 < 0.05的基因数: {sum(corr_df['pval'] < 0.05)}")
    print(f"相关性FDR校正后 < 0.05的基因数: {sum(corr_df['padj'] < 0.05)}")
    print(f"|rho| > 0.3的基因数: {sum(abs(corr_df['rho']) > 0.3)}")
    print(f"|rho| > 0.5的基因数: {sum(abs(corr_df['rho']) > 0.5)}")

def main():
    """
    主分析流程
    """
    parser = argparse.ArgumentParser(description="BRM基因RNA-seq数据分析脚本")
    parser.add_argument("--input_file", "-i", 
                        default="data/processed/working.csv", 
                        help="输入数据文件的路径 (相对于项目根目录)")
    parser.add_argument("--output_dir", "-o", default=None, help="输出目录路径。如果未指定，将使用 'results/<脚本名>'")
    
    args = parser.parse_args()

    # Determine paths based on project root
    input_path = PROJECT_ROOT / args.input_file
    
    if args.output_dir:
        output_path = PROJECT_ROOT / args.output_dir
    else:
        script_name = Path(__file__).stem
        output_path = PROJECT_ROOT / "results" / script_name
    
    # Ensure output directory exists and set up logging
    output_path.mkdir(parents=True, exist_ok=True)
    log_file_path = output_path / "analysis_log.txt"

    original_stdout = sys.stdout
    log_file = open(log_file_path, 'w', encoding='utf-8')
    sys.stdout = Tee(original_stdout, log_file)

    try:
        print("=== BRM基因RNA-seq数据分析 ===\n")
        print(f"输入文件: {input_path}")
        print(f"输出目录: {output_path}")

        if not input_path.exists():
            print(f"错误: 输入文件未找到 at '{input_path}'")
            return

        # 1. 数据加载和预处理
        expr, brm_gene = load_and_preprocess_data(str(input_path))
        
        if expr is None:
            print("数据加载失败，退出分析")
            return
        
        # 2. 计算BRM表达值，进行可视化，并创建分组
        brm_expr = np.log2(expr.loc[brm_gene] + 1)
        
        plot_brm_distribution(brm_expr, brm_gene, output_path)
        
        groups = create_brm_groups(brm_expr, brm_gene)
        
        # 3. 筛选高低表达组的样本进行下游分析
        samples_to_keep = groups[groups.isin(["High", "Low"])].index
        expr_for_analysis = expr[samples_to_keep]
        groups_for_analysis = groups.loc[samples_to_keep]
        print(f"\n为下游分析筛选样本，仅保留高低表达组，共 {len(samples_to_keep)} 个样本。")
        
        # 4. 差异表达分析
        deg_df, significant_genes = differential_expression_analysis(expr_for_analysis, groups_for_analysis, brm_gene)
        
        # 5. 相关性分析
        corr_df, significant_corr = correlation_analysis(expr_for_analysis, brm_gene)
        
        # 6. 保存结果
        save_results(deg_df, significant_genes, corr_df, significant_corr, output_path)
        
        # 7. 输出分析总结
        print_summary(brm_gene, expr_for_analysis, significant_genes, significant_corr, deg_df, corr_df)

    finally:
        # Restore stdout and close the log file
        sys.stdout = original_stdout
        log_file.close()
        print(f"\n分析日志已保存到: {log_file_path}")

if __name__ == "__main__":
    main() 