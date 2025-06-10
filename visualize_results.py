"""
BRM基因RNA-seq数据分析结果可视化脚本
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# 设置专业出版级图形样式
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 16
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['xtick.minor.width'] = 0.8
plt.rcParams['ytick.minor.width'] = 0.8
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linewidth'] = 0.8

# 设置专业配色方案
colors_pub = {
    'red': '#E74C3C',      # 专业红色
    'blue': '#3498DB',     # 专业蓝色  
    'green': '#2ECC71',    # 专业绿色
    'orange': '#F39C12',   # 专业橙色
    'purple': '#9B59B6',   # 专业紫色
    'gray': '#95A5A6',     # 专业灰色
    'dark_red': '#C0392B',
    'dark_blue': '#2980B9',
    'light_blue': '#AED6F1',
    'light_red': '#F1948A'
}

sns.set_style("whitegrid", {
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.2,
    "grid.linewidth": 0.8,
    "font.family": ["Arial"]
})

def load_data():
    """
    Load analysis results data
    """
    print("Loading analysis results data...")
    
    # Load differential expression results
    deg_df = pd.read_csv('differential_expression_results.csv')
    significant_genes = pd.read_csv('significant_DEGs.csv')
    
    # Load correlation results
    corr_df = pd.read_csv('correlation_analysis_results.csv')
    significant_corr = pd.read_csv('significant_correlations.csv')
    
    # Load original expression data
    expr = pd.read_csv("illuminahiseq_rnaseqv2-RSEM_genes_normalized (MD5).csv", index_col=0)
    
    # Data cleaning
    if 'Hybridization REF' in expr.index:
        expr = expr.drop('Hybridization REF')
    
    non_numeric_rows = []
    for idx in expr.index:
        if idx in ['gene_id', 'normalized_count', 'Hybridization REF']:
            non_numeric_rows.append(idx)
    
    if non_numeric_rows:
        expr = expr.drop(non_numeric_rows)
    
    expr = expr.apply(pd.to_numeric, errors='coerce')
    expr = expr.dropna()
    
    # Filter low expression genes
    expr_filtered = expr[(expr > 1).sum(axis=1) >= 0.1 * expr.shape[1]]
    
    print(f"Number of DEGs: {len(significant_genes)}")
    print(f"Number of significantly correlated genes: {len(significant_corr)}")
    print(f"Expression data dimensions: {expr_filtered.shape}")
    
    return deg_df, significant_genes, corr_df, significant_corr, expr_filtered

def create_volcano_plot(deg_df):
    """
    Create volcano plot
    """
    print("Creating volcano plot...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Calculate -log10(padj)
    deg_df['-log10_padj'] = -np.log10(deg_df['padj'] + 1e-300)
    
    # Create color mapping based on significance and fold change
    def get_point_color(row):
        if row['padj'] < 0.05:
            if row['log2FC'] > 1:
                return colors_pub['red']
            elif row['log2FC'] < -1:
                return colors_pub['blue']
            else:
                return colors_pub['orange']
        else:
            return colors_pub['gray']
    
    point_colors = [get_point_color(row) for _, row in deg_df.iterrows()]
    
    # Create scatter plot with different sizes for significance
    sizes = [12 if row['padj'] < 0.001 else 8 if row['padj'] < 0.05 else 4 
             for _, row in deg_df.iterrows()]
    
    scatter = ax.scatter(deg_df['log2FC'], deg_df['-log10_padj'], 
                        c=point_colors, alpha=0.7, s=sizes, 
                        edgecolors='none', rasterized=True)
    
    # Add threshold lines with better styling
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.7, linewidth=1.5, zorder=0)
    ax.axvline(x=1, color='black', linestyle='--', alpha=0.7, linewidth=1.5, zorder=0)
    ax.axvline(x=-1, color='black', linestyle='--', alpha=0.7, linewidth=1.5, zorder=0)
    
    # Annotate top genes with better positioning
    top_genes = deg_df.nlargest(8, '-log10_padj')
    for _, gene in top_genes.iterrows():
        if abs(gene['log2FC']) > 0.8 and gene['-log10_padj'] > 5:
            ax.annotate(gene['gene'].split('|')[0], 
                       (gene['log2FC'], gene['-log10_padj']),
                       xytext=(8, 8), textcoords='offset points',
                       fontsize=9, alpha=0.9, fontweight='bold',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor='white', 
                                edgecolor='gray', alpha=0.8),
                       arrowprops=dict(arrowstyle='->', color='gray', alpha=0.6))
    
    # Set labels and title
    ax.set_xlabel('log₂(Fold Change)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log₁₀(adjusted p-value)', fontsize=12, fontweight='bold')
    ax.set_title('Volcano Plot: SMARCA2 High vs Low Expression Groups', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add statistics box with better formatting
    n_up = len(deg_df[(deg_df['padj'] < 0.05) & (deg_df['log2FC'] > 0)])
    n_down = len(deg_df[(deg_df['padj'] < 0.05) & (deg_df['log2FC'] < 0)])
    n_total = len(deg_df[deg_df['padj'] < 0.05])
    
    stats_text = f'Significant DEGs: {n_total}\nUpregulated: {n_up}\nDownregulated: {n_down}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', fontweight='bold',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                    edgecolor='black', alpha=0.9))
    
    # Create professional legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors_pub['red'], label=f'Significantly Up (|log₂FC|>1, n={len(deg_df[(deg_df["padj"] < 0.05) & (deg_df["log2FC"] > 1)])})'),
        Patch(facecolor=colors_pub['blue'], label=f'Significantly Down (|log₂FC|>1, n={len(deg_df[(deg_df["padj"] < 0.05) & (deg_df["log2FC"] < -1)])})'),
        Patch(facecolor=colors_pub['orange'], label=f'Significant (|log₂FC|<1, n={len(deg_df[(deg_df["padj"] < 0.05) & (abs(deg_df["log2FC"]) < 1)])})'),
        Patch(facecolor=colors_pub['gray'], label=f'Not significant (n={len(deg_df[deg_df["padj"] >= 0.05])})')
    ]
    
    legend = ax.legend(handles=legend_elements, loc='upper right', 
                      frameon=True, fancybox=True, shadow=True,
                      framealpha=0.9, fontsize=9)
    legend.get_frame().set_linewidth(1.2)
    
    # Set axis limits with padding
    x_margin = max(abs(deg_df['log2FC'])) * 0.1
    y_margin = max(deg_df['-log10_padj']) * 0.05
    ax.set_xlim(deg_df['log2FC'].min() - x_margin, deg_df['log2FC'].max() + x_margin)
    ax.set_ylim(-0.5, deg_df['-log10_padj'].max() + y_margin)
    
    # Improve grid
    ax.grid(True, alpha=0.3, linewidth=0.8)
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    plt.savefig('fig/volcano_plot.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig('fig/volcano_plot.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig('fig/volcano_plot.eps', bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

def create_expression_distribution(expr_filtered):
    """
    Create SMARCA2 expression distribution plot
    """
    print("Creating SMARCA2 expression distribution...")
    
    # Find SMARCA2 gene
    smarca2_gene = None
    for gene in expr_filtered.index:
        if 'SMARCA2' in str(gene).upper():
            smarca2_gene = gene
            break
    
    if smarca2_gene is None:
        print("SMARCA2 gene not found")
        return
    
    # Get expression data
    smarca2_expr = np.log2(expr_filtered.loc[smarca2_gene] + 1)
    smarca2_median = np.median(smarca2_expr)
    smarca2_mean = np.mean(smarca2_expr)
    groups = np.where(smarca2_expr > smarca2_median, "High", "Low")
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Histogram with better styling
    n, bins, patches = axes[0].hist(smarca2_expr, bins=25, alpha=0.8, 
                                   color=colors_pub['blue'], edgecolor='white', linewidth=1.2)
    
    # Add median and mean lines
    axes[0].axvline(smarca2_median, color=colors_pub['red'], linestyle='--', linewidth=2.5, 
                   label=f'Median: {smarca2_median:.2f}', alpha=0.9)
    axes[0].axvline(smarca2_mean, color=colors_pub['orange'], linestyle=':', linewidth=2.5, 
                   label=f'Mean: {smarca2_mean:.2f}', alpha=0.9)
    
    axes[0].set_xlabel('SMARCA2 Expression (log₂(TPM+1))', fontweight='bold')
    axes[0].set_ylabel('Number of Samples', fontweight='bold')
    axes[0].set_title('SMARCA2 Expression Distribution', fontweight='bold', pad=15)
    
    legend1 = axes[0].legend(frameon=True, fancybox=True, shadow=True, framealpha=0.9)
    legend1.get_frame().set_linewidth(1.2)
    axes[0].grid(True, alpha=0.3, linewidth=0.8)
    axes[0].set_axisbelow(True)
    
    # Box plot with better styling
    box_data = [smarca2_expr[groups == 'Low'], smarca2_expr[groups == 'High']]
    bp = axes[1].boxplot(box_data, labels=['Low Expression\nGroup', 'High Expression\nGroup'], 
                        patch_artist=True, widths=0.6)
    
    # Customize box plot colors
    bp['boxes'][0].set_facecolor(colors_pub['light_blue'])
    bp['boxes'][1].set_facecolor(colors_pub['light_red'])
    bp['boxes'][0].set_edgecolor(colors_pub['dark_blue'])
    bp['boxes'][1].set_edgecolor(colors_pub['dark_red'])
    bp['boxes'][0].set_linewidth(1.5)
    bp['boxes'][1].set_linewidth(1.5)
    
    # Customize whiskers, caps and medians
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black', linewidth=1.5)
    
    axes[1].set_ylabel('SMARCA2 Expression (log₂(TPM+1))', fontweight='bold')
    axes[1].set_title('SMARCA2 Expression by Groups', fontweight='bold', pad=15)
    axes[1].grid(True, alpha=0.3, linewidth=0.8)
    axes[1].set_axisbelow(True)
    
    # Add sample size annotations
    n_low = len(box_data[0])
    n_high = len(box_data[1])
    axes[1].text(1, axes[1].get_ylim()[1]*0.95, f'n={n_low}', ha='center', fontweight='bold')
    axes[1].text(2, axes[1].get_ylim()[1]*0.95, f'n={n_high}', ha='center', fontweight='bold')
    
    # Violin plot instead of simple density plot
    violin_parts = axes[2].violinplot([smarca2_expr[groups == 'Low'], smarca2_expr[groups == 'High']], 
                                     positions=[1, 2], widths=0.8, showmeans=True, showmedians=True)
    
    # Customize violin plot
    for i, pc in enumerate(violin_parts['bodies']):
        if i == 0:
            pc.set_facecolor(colors_pub['light_blue'])
            pc.set_edgecolor(colors_pub['dark_blue'])
        else:
            pc.set_facecolor(colors_pub['light_red'])
            pc.set_edgecolor(colors_pub['dark_red'])
        pc.set_alpha(0.8)
        pc.set_linewidth(1.5)
    
    # Customize other violin elements
    for element in ['cbars', 'cmins', 'cmaxes', 'cmedians', 'cmeans']:
        if element in violin_parts:
            violin_parts[element].set_color('black')
            violin_parts[element].set_linewidth(1.5)
    
    axes[2].set_xticks([1, 2])
    axes[2].set_xticklabels(['Low Expression\nGroup', 'High Expression\nGroup'])
    axes[2].set_ylabel('SMARCA2 Expression (log₂(TPM+1))', fontweight='bold')
    axes[2].set_title('SMARCA2 Expression Density Distribution', fontweight='bold', pad=15)
    axes[2].grid(True, alpha=0.3, linewidth=0.8)
    axes[2].set_axisbelow(True)
    
    # Add overall statistics
    from scipy import stats
    statistic, p_value = stats.mannwhitneyu(box_data[0], box_data[1], alternative='two-sided')
    
    fig.suptitle(f'SMARCA2 Expression Analysis (Mann-Whitney U test p-value: {p_value:.2e})', 
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig('fig/smarca2_expression_distribution.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig('fig/smarca2_expression_distribution.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig('fig/smarca2_expression_distribution.eps', bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

def create_deg_heatmap(expr_filtered, significant_genes, top_n=50):
    """
    Create differentially expressed genes heatmap
    """
    print(f"Creating Top {top_n} DEGs heatmap...")
    
    if len(significant_genes) == 0:
        print("No significant DEGs found")
        return
    
    # 选择top基因
    top_genes_list = significant_genes.head(top_n)["gene"].tolist()
    
    # 检查基因是否在表达数据中
    available_genes = [g for g in top_genes_list if g in expr_filtered.index]
    if len(available_genes) < 10:
        print(f"Only {len(available_genes)} genes available in expression data, skipping heatmap")
        return
    
    # 准备数据
    heatmap_data = np.log2(expr_filtered.loc[available_genes] + 1)
    
    # 创建分组信息
    smarca2_gene = None
    for gene in expr_filtered.index:
        if 'SMARCA2' in str(gene).upper():
            smarca2_gene = gene
            break
    
    if smarca2_gene is None:
        print("SMARCA2 gene not found for grouping")
        return
    
    smarca2_expr = np.log2(expr_filtered.loc[smarca2_gene] + 1)
    smarca2_median = np.median(smarca2_expr)
    groups = np.where(smarca2_expr > smarca2_median, "High", "Low")
    groups_series = pd.Series(groups, index=expr_filtered.columns)
    
    # 准备列颜色
    col_colors = groups_series.map({"High": "red", "Low": "blue"})
    
    # 创建基因标签（只保留基因符号）
    gene_labels = [g.split('|')[0] if '|' in g else g for g in available_genes]
    
    plt.figure(figsize=(15, 12))
    
    try:
        # 绘制聚类热图
        g = sns.clustermap(heatmap_data, 
                          col_colors=col_colors,
                          figsize=(15, 12),
                          cmap='RdBu_r',
                          center=0,
                          xticklabels=False,
                          yticklabels=gene_labels,
                          cbar_kws={'label': 'log2(TPM+1)'},
                          row_cluster=True,
                          col_cluster=True)
        
        # 添加标题
        g.fig.suptitle(f'Top {len(available_genes)} Differentially Expressed Genes Heatmap', 
                      fontsize=16, y=0.95)
        
        g.savefig('fig/deg_heatmap.png', dpi=600, bbox_inches='tight', facecolor='white')
        g.savefig('fig/deg_heatmap.pdf', bbox_inches='tight', facecolor='white')
        g.savefig('fig/deg_heatmap.eps', bbox_inches='tight', facecolor='white')
        plt.show()
        
    except Exception as e:
        print(f"Heatmap creation failed: {e}")
        # Create simple heatmap as backup
        plt.figure(figsize=(15, 10))
        sns.heatmap(heatmap_data, 
                   cmap='RdBu_r', center=0,
                   xticklabels=False,
                   yticklabels=gene_labels,
                   cbar_kws={'label': 'log2(TPM+1)'})
        plt.title(f'Top {len(available_genes)} Differentially Expressed Genes')
        plt.tight_layout()
        plt.savefig('fig/deg_heatmap_simple.png', dpi=600, bbox_inches='tight', facecolor='white')
        plt.show()
        plt.close()

def create_correlation_plots(corr_df, significant_corr):
    """
    Create correlation analysis plots
    """
    print("Creating correlation analysis plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Correlation coefficient distribution histogram
    axes[0,0].hist(corr_df['rho'], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    axes[0,0].axvline(0, color='red', linestyle='--', alpha=0.7)
    axes[0,0].set_xlabel('Spearman Correlation Coefficient')
    axes[0,0].set_ylabel('Number of Genes')
    axes[0,0].set_title('Correlation Coefficient Distribution with SMARCA2')
    axes[0,0].grid(True, alpha=0.3)
    
    # 2. Significantly correlated genes scatter plot
    sig_pos = significant_corr[significant_corr['rho'] > 0]
    sig_neg = significant_corr[significant_corr['rho'] < 0]
    
    axes[0,1].scatter(range(len(sig_pos)), sig_pos['rho'], 
                     color='red', alpha=0.6, s=20, label=f'Positive correlation (n={len(sig_pos)})')
    axes[0,1].scatter(range(len(sig_neg)), sig_neg['rho'], 
                     color='blue', alpha=0.6, s=20, label=f'Negative correlation (n={len(sig_neg)})')
    axes[0,1].axhline(0.3, color='red', linestyle='--', alpha=0.5)
    axes[0,1].axhline(-0.3, color='blue', linestyle='--', alpha=0.5)
    axes[0,1].set_xlabel('Gene Rank')
    axes[0,1].set_ylabel('Correlation Coefficient')
    axes[0,1].set_title('Significantly Correlated Genes')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # 3. P-value distribution
    axes[1,0].hist(-np.log10(corr_df['pval'] + 1e-300), bins=50, 
                  alpha=0.7, color='green', edgecolor='black')
    axes[1,0].axvline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7, 
                     label='p=0.05')
    axes[1,0].set_xlabel('-log10(p-value)')
    axes[1,0].set_ylabel('Number of Genes')
    axes[1,0].set_title('Correlation P-value Distribution')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)
    
    # 4. Correlation coefficient vs P-value
    colors = ['red' if abs(rho) > 0.3 and padj < 0.05 else 'gray' 
              for rho, padj in zip(corr_df['rho'], corr_df['padj'])]
    
    axes[1,1].scatter(corr_df['rho'], -np.log10(corr_df['padj'] + 1e-300),
                     c=colors, alpha=0.6, s=8)
    axes[1,1].axvline(0.3, color='red', linestyle='--', alpha=0.5)
    axes[1,1].axvline(-0.3, color='red', linestyle='--', alpha=0.5)
    axes[1,1].axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5)
    axes[1,1].set_xlabel('Correlation Coefficient')
    axes[1,1].set_ylabel('-log10(adjusted p-value)')
    axes[1,1].set_title('Correlation Coefficient vs Significance')
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('fig/correlation_analysis.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig('fig/correlation_analysis.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig('fig/correlation_analysis.eps', bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

def create_top_genes_barplot(significant_genes, top_n=15):
    """
    Create top DEGs bar plot
    """
    print(f"Creating Top {top_n} DEGs bar plot...")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 10))
    
    # Top upregulated genes
    top_up = significant_genes[significant_genes['log2FC'] > 0].head(top_n)
    if len(top_up) > 0:
        gene_names = [g.split('|')[0] if '|' in g else g for g in top_up['gene']]
        y_pos = np.arange(len(gene_names))
        
        # Create gradient colors based on fold change
        fc_values = top_up['log2FC'].values
        colors_up = [plt.cm.Reds(0.5 + 0.5 * (fc / max(fc_values))) for fc in fc_values]
        
        bars1 = axes[0].barh(y_pos, top_up['log2FC'], color=colors_up, 
                            alpha=0.8, edgecolor='black', linewidth=0.8)
        
        axes[0].set_yticks(y_pos)
        axes[0].set_yticklabels(gene_names, fontsize=10, fontweight='bold')
        axes[0].set_xlabel('log₂(Fold Change)', fontsize=12, fontweight='bold')
        axes[0].set_title(f'Top {len(top_up)} Upregulated Genes', 
                         fontsize=14, fontweight='bold', pad=20)
        axes[0].grid(True, alpha=0.3, linewidth=0.8, axis='x')
        axes[0].set_axisbelow(True)
        
        # Add fold change values on bars
        for i, (bar, fc, pval) in enumerate(zip(bars1, top_up['log2FC'], top_up['padj'])):
            width = bar.get_width()
            axes[0].text(width + 0.05, bar.get_y() + bar.get_height()/2, 
                        f'{fc:.2f}', va='center', fontsize=9, fontweight='bold')
        
        # Set x-axis limits with padding
        axes[0].set_xlim(0, max(top_up['log2FC']) * 1.3)
        
        # Add significance indicators
        for i, pval in enumerate(top_up['padj']):
            if pval < 0.001:
                sig_text = '***'
            elif pval < 0.01:
                sig_text = '**'
            elif pval < 0.05:
                sig_text = '*'
            else:
                sig_text = 'ns'
            
            axes[0].text(-0.1, i, sig_text, va='center', ha='right', 
                        fontsize=10, fontweight='bold', color=colors_pub['red'])
    
    # Top downregulated genes
    top_down = significant_genes[significant_genes['log2FC'] < 0].head(top_n)
    if len(top_down) > 0:
        gene_names = [g.split('|')[0] if '|' in g else g for g in top_down['gene']]
        y_pos = np.arange(len(gene_names))
        
        # Create gradient colors based on fold change (reversed for downregulation)
        fc_values = top_down['log2FC'].values
        colors_down = [plt.cm.Blues(0.5 + 0.5 * (abs(fc) / max(abs(fc_values)))) for fc in fc_values]
        
        bars2 = axes[1].barh(y_pos, top_down['log2FC'], color=colors_down, 
                            alpha=0.8, edgecolor='black', linewidth=0.8)
        
        axes[1].set_yticks(y_pos)
        axes[1].set_yticklabels(gene_names, fontsize=10, fontweight='bold')
        axes[1].set_xlabel('log₂(Fold Change)', fontsize=12, fontweight='bold')
        axes[1].set_title(f'Top {len(top_down)} Downregulated Genes', 
                         fontsize=14, fontweight='bold', pad=20)
        axes[1].grid(True, alpha=0.3, linewidth=0.8, axis='x')
        axes[1].set_axisbelow(True)
        
        # Add fold change values on bars
        for i, (bar, fc, pval) in enumerate(zip(bars2, top_down['log2FC'], top_down['padj'])):
            width = bar.get_width()
            axes[1].text(width - 0.05, bar.get_y() + bar.get_height()/2, 
                        f'{fc:.2f}', va='center', ha='right', fontsize=9, fontweight='bold')
        
        # Set x-axis limits with padding
        axes[1].set_xlim(min(top_down['log2FC']) * 1.3, 0)
        
        # Add significance indicators
        for i, pval in enumerate(top_down['padj']):
            if pval < 0.001:
                sig_text = '***'
            elif pval < 0.01:
                sig_text = '**'
            elif pval < 0.05:
                sig_text = '*'
            else:
                sig_text = 'ns'
            
            axes[1].text(0.1, i, sig_text, va='center', ha='left', 
                        fontsize=10, fontweight='bold', color=colors_pub['blue'])
    
    # Add significance legend
    from matplotlib.patches import Patch
    sig_legend = [
        Patch(facecolor='none', label='Significance:'),
        Patch(facecolor='none', label='*** p < 0.001'),
        Patch(facecolor='none', label='** p < 0.01'),
        Patch(facecolor='none', label='* p < 0.05'),
        Patch(facecolor='none', label='ns p ≥ 0.05')
    ]
    
    fig.legend(handles=sig_legend, loc='center', bbox_to_anchor=(0.5, 0.02), 
              ncol=5, frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig('fig/top_genes_barplot.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig('fig/top_genes_barplot.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig('fig/top_genes_barplot.eps', bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

def create_summary_stats_plot(deg_df, significant_genes, significant_corr):
    """
    Create analysis results summary plot
    """
    print("Creating analysis results summary...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. DEG statistics
    up_genes = len(significant_genes[significant_genes['log2FC'] > 0])
    down_genes = len(significant_genes[significant_genes['log2FC'] < 0])
    
    axes[0,0].pie([up_genes, down_genes], 
                 labels=[f'Upregulated ({up_genes})', f'Downregulated ({down_genes})'],
                 colors=['red', 'blue'], autopct='%1.1f%%')
    axes[0,0].set_title('Significant DEGs Distribution')
    
    # 2. Effect size distribution
    fc_bins = [0, 0.5, 1, 1.5, 2, float('inf')]
    fc_labels = ['0-0.5', '0.5-1', '1-1.5', '1.5-2', '>2']
    fc_counts = []
    
    for i in range(len(fc_bins)-1):
        count = len(significant_genes[(abs(significant_genes['log2FC']) >= fc_bins[i]) & 
                                    (abs(significant_genes['log2FC']) < fc_bins[i+1])])
        fc_counts.append(count)
    
    axes[0,1].bar(fc_labels, fc_counts, color='skyblue', alpha=0.7)
    axes[0,1].set_xlabel('|log2(Fold Change)|')
    axes[0,1].set_ylabel('Number of Genes')
    axes[0,1].set_title('Effect Size Distribution')
    axes[0,1].tick_params(axis='x', rotation=45)
    
    # 3. Significance level distribution
    p_bins = [0, 0.001, 0.01, 0.05, 1]
    p_labels = ['<0.001', '0.001-0.01', '0.01-0.05', '0.05-1']
    p_counts = []
    
    for i in range(len(p_bins)-1):
        count = len(deg_df[(deg_df['padj'] >= p_bins[i]) & 
                          (deg_df['padj'] < p_bins[i+1])])
        p_counts.append(count)
    
    axes[1,0].bar(p_labels, p_counts, color='lightgreen', alpha=0.7)
    axes[1,0].set_xlabel('Adjusted P-value')
    axes[1,0].set_ylabel('Number of Genes')
    axes[1,0].set_title('Significance Level Distribution')
    axes[1,0].tick_params(axis='x', rotation=45)
    
    # 4. Correlation statistics
    pos_corr = len(significant_corr[significant_corr['rho'] > 0])
    neg_corr = len(significant_corr[significant_corr['rho'] < 0])
    
    axes[1,1].pie([pos_corr, neg_corr], 
                 labels=[f'Positive correlation ({pos_corr})', f'Negative correlation ({neg_corr})'],
                 colors=['orange', 'purple'], autopct='%1.1f%%')
    axes[1,1].set_title('Significant Correlated Genes Distribution')
    
    plt.tight_layout()
    plt.savefig('fig/summary_statistics.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig('fig/summary_statistics.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig('fig/summary_statistics.eps', bbox_inches='tight', facecolor='white')
    plt.show()
    plt.close()

def main():
    """
    Main visualization workflow
    """
    print("=== BRM Gene RNA-seq Data Analysis Results Visualization ===\n")
    
    # Load data
    deg_df, significant_genes, corr_df, significant_corr, expr_filtered = load_data()
    
    # Generate various plots
    print("\n1. Creating volcano plot...")
    create_volcano_plot(deg_df)
    
    print("\n2. Creating expression distribution plot...")
    create_expression_distribution(expr_filtered)
    
    print("\n3. Creating DEGs heatmap...")
    create_deg_heatmap(expr_filtered, significant_genes, top_n=50)
    
    print("\n4. Creating correlation analysis plots...")
    create_correlation_plots(corr_df, significant_corr)
    
    print("\n5. Creating top genes bar plot...")
    create_top_genes_barplot(significant_genes, top_n=15)
    
    print("\n6. Creating summary statistics plot...")
    create_summary_stats_plot(deg_df, significant_genes, significant_corr)
    
    print("\n=== Visualization Complete ===")
    print("All plots have been saved to fig/ directory in both PNG and PDF formats")
    print("Generated plots include:")
    print("- volcano_plot: Volcano plot")
    print("- smarca2_expression_distribution: SMARCA2 expression distribution")
    print("- deg_heatmap: Differentially expressed genes heatmap")
    print("- correlation_analysis: Correlation analysis")
    print("- top_genes_barplot: Top genes bar plot")
    print("- summary_statistics: Summary statistics")

if __name__ == "__main__":
    main() 