"""
BRM基因RNA-seq数据分析结果可视化脚本
"""
import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings('ignore')

# Define project root
PROJECT_ROOT = Path(__file__).resolve().parents[1]

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
    'red': '#E74C3C',
    'blue': '#3498DB',
    'green': '#2ECC71',
    'orange': '#F39C12',
    'purple': '#9B59B6',
    'gray': '#95A5A6',
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


def load_data(input_dir: Path, expression_file: Path, target_gene: str):
    """
    Load analysis results data.
    """
    print("Loading analysis results data...")

    try:
        # Use fixed filenames as found in the output directory
        deg_path = input_dir / "differential_expression_results.csv"
        sig_genes_path = input_dir / "significant_DEGs.csv"
        corr_path = input_dir / "correlation_analysis_results.csv"
        sig_corr_path = input_dir / "significant_correlations.csv"

        print(f" - Loading DEGs from: {deg_path.relative_to(PROJECT_ROOT)}")
        deg_df = pd.read_csv(deg_path)

        print(f" - Loading significant DEGs from: {sig_genes_path.relative_to(PROJECT_ROOT)}")
        significant_genes = pd.read_csv(sig_genes_path)

        print(f" - Loading correlations from: {corr_path.relative_to(PROJECT_ROOT)}")
        corr_df = pd.read_csv(corr_path)

        print(f" - Loading significant correlations from: {sig_corr_path.relative_to(PROJECT_ROOT)}")
        significant_corr = pd.read_csv(sig_corr_path)

        print(f" - Loading expression data from: {expression_file.relative_to(PROJECT_ROOT)}")
        expr = pd.read_csv(expression_file, index_col=0)

    except FileNotFoundError as e:
        print(f"\nError loading data: {e}")
        print(
            f"Please ensure that the analysis for '{target_gene}' has been run and results are in '{input_dir.relative_to(PROJECT_ROOT)}'")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), None

    # Filter low expression genes, a common practice for visualization clarity
    expr_filtered = expr[(expr > 1).sum(axis=1) >= 0.1 * expr.shape[1]]

    # Find full gene name (e.g., SMARCA2|6595) in the expression data
    target_gene_full_name = None
    for gene in expr_filtered.index:
        if str(gene).upper().startswith(target_gene.upper()):
            target_gene_full_name = gene
            print(f"Found target gene in expression data: {target_gene_full_name}")
            break

    if target_gene_full_name is None:
        print(f"Error: Target gene '{target_gene}' not found in the expression data index.")
    else:
        # Exclude the target gene itself from the results
        print(f"Excluding target gene '{target_gene_full_name}' from analysis results...")
        if 'gene' in deg_df.columns:
            deg_df = deg_df[deg_df['gene'] != target_gene_full_name]
        if 'gene' in significant_genes.columns:
            significant_genes = significant_genes[significant_genes['gene'] != target_gene_full_name]
        if 'gene' in corr_df.columns:
            corr_df = corr_df[corr_df['gene'] != target_gene_full_name]
        if 'gene' in significant_corr.columns:
            significant_corr = significant_corr[significant_corr['gene'] != target_gene_full_name]
        print("  - Target gene excluded from DEG and correlation tables.")

    print(f"\nSuccessfully loaded data:")
    print(f"  - Total DEGs analyzed: {len(deg_df)}")
    print(f"  - Significant DEGs found: {len(significant_genes)}")
    print(f"  - Total correlations analyzed: {len(corr_df)}")
    print(f"  - Significant correlations found: {len(significant_corr)}")
    print(f"  - Expression data dimensions after filtering: {expr_filtered.shape}")

    return deg_df, significant_genes, corr_df, significant_corr, expr_filtered, target_gene_full_name


def create_volcano_plot(deg_df: pd.DataFrame, output_dir: Path, target_gene: str):
    """
    Create volcano plot, styled for publication.
    """
    print("Creating volcano plot...")

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.despine(ax=ax) # More professional look

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

    ax.scatter(deg_df['log2FC'], deg_df['-log10_padj'],
               c=point_colors, alpha=0.7, s=sizes,
               edgecolors='none', rasterized=True)

    # Add threshold lines with better styling
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.7, linewidth=1.5, zorder=0)
    ax.axvline(x=1, color='black', linestyle='--', alpha=0.7, linewidth=1.5, zorder=0)
    ax.axvline(x=-1, color='black', linestyle='--', alpha=0.7, linewidth=1.5, zorder=0)

    # Annotate top genes (top 5 by padj) with predefined offsets to avoid overlap
    top_genes = deg_df.nsmallest(5, 'padj')
    
    # Predefined offsets for labels
    offsets = [(20, 20), (-20, -20), (20, -20), (-20, 20), (30, 0)]

    for i, (_, gene) in enumerate(top_genes.iterrows()):
        ax.annotate(gene['gene'].split('|')[0],
                    (gene['log2FC'], gene['-log10_padj']),
                    xytext=offsets[i % len(offsets)], textcoords='offset points',
                    fontsize=9, alpha=0.9, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='white',
                              edgecolor='gray', alpha=0.8),
                    arrowprops=dict(arrowstyle='->', color='gray', alpha=0.6))

    # Set labels and title
    ax.set_xlabel('log2(Fold Change)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(adjusted p-value)', fontsize=12, fontweight='bold')
    ax.set_title(f'Volcano Plot: {target_gene} High vs Low Expression Groups',
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

    # Create professional, simplified legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors_pub['red'],
              label=f'Upregulated (|log2FC|>1)'),
        Patch(facecolor=colors_pub['blue'],
              label=f'Downregulated (|log2FC|>1)'),
        Patch(facecolor=colors_pub['orange'],
              label=f'Significant (p<0.05)'),
        Patch(facecolor=colors_pub['gray'], label='Not significant')
    ]

    legend = ax.legend(handles=legend_elements, loc='upper right',
                       frameon=False, fontsize=9) # Cleaner legend without frame

    # Set axis limits with padding
    x_margin = max(abs(deg_df['log2FC'])) * 0.1
    y_margin = max(deg_df['-log10_padj']) * 0.05
    ax.set_xlim(deg_df['log2FC'].min() - x_margin, deg_df['log2FC'].max() + x_margin)
    ax.set_ylim(-0.5, deg_df['-log10_padj'].max() + y_margin)

    # Improve grid
    ax.grid(True, alpha=0.3, linewidth=0.8)
    ax.set_axisbelow(True)

    plt.tight_layout()
    output_path = output_dir / 'volcano_plot'
    print(f"  Saving volcano plot to {output_path}...")
    plt.savefig(f'{output_path}.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.eps', bbox_inches='tight', facecolor='white')
    # plt.show()
    plt.close()


def create_expression_distribution(expr_filtered: pd.DataFrame, target_gene_full_name: str, output_dir: Path,
                                   target_gene: str):
    """
    Create target gene expression distribution plot, styled for publication.
    This combines a histogram and a violin+box plot with integrated statistical annotation.
    """
    print(f"Creating {target_gene} expression distribution...")

    if target_gene_full_name is None or target_gene_full_name not in expr_filtered.index:
        print(
            f"Target gene '{target_gene}' ({target_gene_full_name}) not found. Skipping distribution plot.")
        return

    # --- Data Preparation ---
    target_expr = np.log2(expr_filtered.loc[target_gene_full_name] + 1)
    target_median = np.median(target_expr)
    target_mean = np.mean(target_expr)
    groups = pd.Series(np.where(target_expr > target_median, "High", "Low"), index=target_expr.index)
    
    low_expr = target_expr[groups == 'Low']
    high_expr = target_expr[groups == 'High']
    
    # --- Statistical Analysis ---
    from scipy import stats
    statistic, p_value = stats.mannwhitneyu(low_expr, high_expr, alternative='two-sided')
    if p_value < 0.001: p_text = '***'
    elif p_value < 0.01: p_text = '**'
    elif p_value < 0.05: p_text = '*'
    else: p_text = f'p = {p_value:.2f}'

    # --- Plotting (2-panel layout) ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [1, 1.2]})
    fig.suptitle(f'{target_gene} Expression Analysis', fontsize=16, fontweight='bold')

    # Panel A: Histogram
    ax1 = axes[0]
    sns.histplot(target_expr, bins=25, ax=ax1, color=colors_pub['blue'], alpha=0.7, edgecolor='white')
    ax1.axvline(target_median, color=colors_pub['red'], linestyle='--', linewidth=2, label=f'Median: {target_median:.2f}')
    ax1.set_title('A. Overall Distribution', loc='left', fontsize=14, fontweight='bold')
    ax1.set_xlabel(f'{target_gene} Expression (log2(TPM+1))', fontweight='bold')
    ax1.set_ylabel('Number of Samples', fontweight='bold')
    ax1.legend()
    sns.despine(ax=ax1)

    # Panel B: Violin + Box Plot with Stats
    ax2 = axes[1]
    plot_data = pd.DataFrame({'Expression': target_expr, 'Group': groups})
    sns.violinplot(x='Group', y='Expression', data=plot_data, order=['Low', 'High'], ax=ax2,
                   palette=[colors_pub['light_blue'], colors_pub['light_red']], inner=None, linewidth=1.5)
    
    # Overlay a stylish boxplot inside the violin
    sns.boxplot(x='Group', y='Expression', data=plot_data, order=['Low', 'High'], ax=ax2,
                width=0.2, boxprops={'zorder': 2, 'facecolor': 'white'},
                whiskerprops={'linewidth': 1.5},
                capprops={'linewidth': 1.5},
                medianprops={'color': 'black', 'linewidth': 1.5})

    ax2.set_title('B. Comparison by Expression Group', loc='left', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Median-based Group', fontweight='bold')
    ax2.set_ylabel('') # Y-axis label is shared
    ax2.set_xticklabels([f'Low (n={len(low_expr)})', f'High (n={len(high_expr)})'])
    sns.despine(ax=ax2)

    # Add statistical annotation bracket
    y_max = target_expr.max()
    bracket_y = y_max * 1.05
    text_y = y_max * 1.1
    ax2.plot([0, 0, 1, 1], [bracket_y, bracket_y, bracket_y, bracket_y], lw=1.5, c='k')
    ax2.text(0.5, text_y, p_text, ha='center', va='bottom', color='k', fontsize=12)
    ax2.set_ylim(bottom=target_expr.min()*0.9, top=y_max * 1.2)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    output_path = output_dir / f'{target_gene.lower()}_expression_distribution'
    print(f"  Saving expression distribution plot to {output_path}...")
    plt.savefig(f'{output_path}.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.eps', bbox_inches='tight', facecolor='white')
    # plt.show()
    plt.close()


def create_deg_heatmap(expr_filtered: pd.DataFrame, significant_genes: pd.DataFrame, target_gene_full_name: str,
                       output_dir: Path, target_gene: str, top_n=50):
    """
    Create differentially expressed genes heatmap, styled for publication.
    This version uses Z-score normalization and a diverging colormap.
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

    # --- Data Preparation: Use Z-score for professional visualization ---
    # Z-score normalization makes expression patterns comparable across genes.
    from scipy.stats import zscore
    heatmap_data = np.log2(expr_filtered.loc[available_genes] + 1)
    heatmap_data_zscored = heatmap_data.apply(zscore, axis=1).fillna(0)


    # Create grouping information
    if target_gene_full_name is None or target_gene_full_name not in expr_filtered.index:
        print(f"Target gene '{target_gene}' not found for grouping. Skipping heatmap.")
        return

    target_expr = np.log2(expr_filtered.loc[target_gene_full_name] + 1)
    target_median = np.median(target_expr)
    groups = np.where(target_expr > target_median, "High", "Low")
    groups_series = pd.Series(groups, index=expr_filtered.columns)

    # 准备列颜色
    col_colors = groups_series.map({"High": "red", "Low": "blue"})

    # 创建基因标签（只保留基因符号）
    gene_labels = [g.split('|')[0] if '|' in g else g for g in available_genes]

    try:
        # Use a diverging colormap and specify Z-score in the color bar
        g = sns.clustermap(heatmap_data_zscored,
                           col_colors=col_colors,
                           figsize=(15, 13),
                           cmap='vlag',  # Diverging colormap is standard for Z-scores
                           center=0,      # Center the colormap at 0 for Z-scores
                           xticklabels=False,
                           yticklabels=gene_labels,
                           row_cluster=False,
                           col_cluster=False,
                           cbar_kws={'label': 'Z-score'}) # Label for the color bar

        # --- 1. Adjust Color Bar ---
        heatmap_pos = g.ax_heatmap.get_position()
        cbar_width = 0.02
        g.cax.set_position([heatmap_pos.x0 - cbar_width - 0.01, heatmap_pos.y0, cbar_width, heatmap_pos.height])
        g.cax.yaxis.set_ticks_position('left')

        # --- 2. Move gene labels to the right ---
        g.ax_heatmap.yaxis.tick_right()
        g.ax_heatmap.yaxis.set_label_position("right")
        g.ax_heatmap.tick_params(axis='y', labelsize=9, rotation=0)

        # --- 3. Adjust Group Legend ---
        from matplotlib.patches import Patch
        handles = [Patch(facecolor='red', label='High Expression'),
                   Patch(facecolor='blue', label='Low Expression')]
        
        # Get the position of the top color bar to place the legend above it
        col_colors_pos = g.ax_col_colors.get_position()
        
        # Place the legend above the heatmap and center it
        legend_x_center = col_colors_pos.x0 + col_colors_pos.width / 2
        legend_y = col_colors_pos.y1 + 0.01  # Leave a small gap above the color bar
        
        g.fig.legend(handles=handles, title=f'{target_gene} Expression Group',
                     bbox_to_anchor=(legend_x_center, legend_y),
                     loc='lower center', ncol=2, frameon=False, fontsize=10, title_fontsize=11)

        # --- 4. Adjust Main Title ---
        g.fig.suptitle(f'Top {len(available_genes)} Differentially Expressed Genes Heatmap',
                       fontsize=16, fontweight='bold', y=legend_y + 0.05)


        output_path = output_dir / 'deg_heatmap'
        print(f"  Saving DEG heatmap to {output_path}...")
        g.savefig(f'{output_path}.png', dpi=600, facecolor='white')
        g.savefig(f'{output_path}.pdf', facecolor='white')
        g.savefig(f'{output_path}.eps', facecolor='white')
        # plt.show() # Remove to prevent GUI popup
        plt.close(g.fig)

    except Exception as e:
        print(f"Heatmap creation failed: {e}")
        # Create simple heatmap as backup
        plt.figure(figsize=(15, 10))
        sns.heatmap(heatmap_data, # Fallback to non-zscored data
                    cmap='Reds',
                    xticklabels=False,
                    yticklabels=gene_labels,
                    cbar_kws={'label': 'log2(TPM+1)'})
        plt.title(f'Top {len(available_genes)} Differentially Expressed Genes (Fallback)')
        plt.tight_layout()
        output_path = output_dir / 'deg_heatmap_simple'
        print(f"  Saving simple DEG heatmap to {output_path}...")
        plt.savefig(f'{output_path}.png', dpi=600, bbox_inches='tight', facecolor='white')
        # plt.show()
        plt.close()


def create_correlation_plots(corr_df: pd.DataFrame, significant_corr: pd.DataFrame, output_dir: Path, target_gene: str):
    """
    Create correlation analysis plots, styled for publication.
    This is a focused 2-panel plot.
    """
    print("Creating correlation analysis plots...")

    # --- Plotting (2-panel layout) ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(f'Correlation Analysis with {target_gene}', fontsize=16, fontweight='bold')

    # --- Panel A: Correlation Volcano Plot ---
    ax1 = axes[0]
    
    # Prepare data for volcano plot
    corr_df_filtered = corr_df.copy()
    corr_df_filtered['-log10_padj'] = -np.log10(corr_df_filtered['padj'] + 1e-300)
    # Cap the value to avoid extreme outliers stretching the axis
    cap = corr_df_filtered['-log10_padj'].quantile(0.999)
    if pd.notna(cap) and cap > 0:
    corr_df_filtered['-log10_padj'] = corr_df_filtered['-log10_padj'].clip(upper=cap)

    # Define colors
    colors = [colors_pub['red'] if abs(rho) > 0.3 and padj < 0.05 else colors_pub['gray']
              for rho, padj in zip(corr_df_filtered['rho'], corr_df_filtered['padj'])]

    ax1.scatter(corr_df_filtered['rho'], corr_df_filtered['-log10_padj'],
                       c=colors, alpha=0.5, s=10, rasterized=True)
    ax1.axvline(0.3, color='black', linestyle='--', alpha=0.7)
    ax1.axvline(-0.3, color='black', linestyle='--', alpha=0.7)
    ax1.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.7)
    
    ax1.set_title('A. Correlation Coefficient vs. Significance', loc='left', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Spearman Correlation (rho)', fontweight='bold')
    ax1.set_ylabel('-log10(Adjusted P-value)', fontweight='bold')
    sns.despine(ax=ax1)
    ax1.grid(True, alpha=0.3)

    # --- Panel B: Top Correlated Genes Bar Plot ---
    ax2 = axes[1]
    top_n = 10
    top_pos = significant_corr[significant_corr['rho'] > 0].nlargest(top_n, 'rho').sort_values('rho', ascending=True)
    top_neg = significant_corr[significant_corr['rho'] < 0].nsmallest(top_n, 'rho').sort_values('rho', ascending=True)
    
    # Combine for plotting
    top_corr_genes = pd.concat([top_neg, top_pos])
    
    if not top_corr_genes.empty:
      gene_names = [g.split('|')[0] for g in top_corr_genes['gene']]
      bar_colors = [colors_pub['blue'] if rho < 0 else colors_pub['red'] for rho in top_corr_genes['rho']]
      bars = ax2.barh(gene_names, top_corr_genes['rho'], color=bar_colors, edgecolor='black', linewidth=0.8)
      ax2.axvline(0, color='black', linewidth=1)
      
      # Add value labels
      for bar in bars:
          width = bar.get_width()
          ha = 'left' if width > 0 else 'right'
          x_pos = width + (0.01 if width > 0 else -0.01)
          ax2.text(x_pos, bar.get_y() + bar.get_height()/2, f'{width:.2f}', va='center', ha=ha, fontsize=9)

    ax2.set_title(f'B. Top {top_n} Correlated Genes', loc='left', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Spearman Correlation (rho)', fontweight='bold')
    sns.despine(ax=ax2)
    ax2.grid(True, alpha=0.3, axis='x')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    output_path = output_dir / 'correlation_analysis'
    print(f"  Saving correlation analysis plots to {output_path}...")
    plt.savefig(f'{output_path}.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.eps', bbox_inches='tight', facecolor='white')
    # plt.show()
    plt.close()


def create_top_genes_barplot(significant_genes: pd.DataFrame, output_dir: Path, top_n=15):
    """
    Create top DEGs bar plot, styled for publication.
    """
    print(f"Creating Top {top_n} DEGs bar plot...")

    fig, axes = plt.subplots(1, 2, figsize=(16, 10), sharey=False) # sharey=False is important
    fig.suptitle(f'Top {top_n} Differentially Expressed Genes', fontsize=16, fontweight='bold')


    # Top upregulated genes
    ax1 = axes[0]
    sns.despine(ax=ax1)
    top_up = significant_genes[significant_genes['log2FC'] > 0].nlargest(top_n, 'log2FC')
    if len(top_up) > 0:
        top_up = top_up.sort_values('log2FC', ascending=True) # Sort for bar plot
        gene_names = [g.split('|')[0] if '|' in g else g for g in top_up['gene']]
        y_pos = np.arange(len(gene_names))

        # Create gradient colors based on fold change
        fc_values = top_up['log2FC'].values
        colors_up = [plt.cm.Reds(0.5 + 0.5 * (fc / max(fc_values))) for fc in fc_values]

        bars1 = ax1.barh(y_pos, top_up['log2FC'], color=colors_up,
                            alpha=0.8, edgecolor='black', linewidth=0.8)

        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(gene_names, fontsize=10, fontweight='bold')
        ax1.set_xlabel('log2(Fold Change)', fontsize=12, fontweight='bold')
        ax1.set_title('A. Upregulated Genes',
                         fontsize=14, fontweight='bold', pad=20, loc='left')
        ax1.grid(True, alpha=0.3, linewidth=0.8, axis='x')
        ax1.set_axisbelow(True)

        # Add fold change values on bars
        for i, (bar, fc, pval) in enumerate(zip(bars1, top_up['log2FC'], top_up['padj'])):
            width = bar.get_width()
            ax1.text(width + 0.05, bar.get_y() + bar.get_height() / 2,
                         f'{fc:.2f}', va='center', fontsize=9, fontweight='bold')

        # Set x-axis limits with padding
        ax1.set_xlim(0, max(top_up['log2FC']) * 1.4)

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

            # Use axis transform to precisely position stars
            transform = ax1.get_yaxis_transform()
            ax1.text(-0.15, i, sig_text, transform=transform, va='center', ha='right',
                         fontsize=10, fontweight='bold', color=colors_pub['red'])

    # Top downregulated genes
    ax2 = axes[1]
    sns.despine(ax=ax2)
    top_down = significant_genes[significant_genes['log2FC'] < 0].nsmallest(top_n, 'log2FC')
    if len(top_down) > 0:
        top_down = top_down.sort_values('log2FC', ascending=True) # Sort for bar plot
        gene_names = [g.split('|')[0] if '|' in g else g for g in top_down['gene']]
        y_pos = np.arange(len(gene_names))

        # Create gradient colors based on fold change (reversed for downregulation)
        fc_values = top_down['log2FC'].values
        colors_down = [plt.cm.Blues(0.5 + 0.5 * (abs(fc) / max(abs(fc_values)))) for fc in fc_values]

        bars2 = ax2.barh(y_pos, top_down['log2FC'], color=colors_down,
                             alpha=0.8, edgecolor='black', linewidth=0.8)

        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(gene_names, fontsize=10, fontweight='bold')
        ax2.set_xlabel('log2(Fold Change)', fontsize=12, fontweight='bold')
        ax2.set_title('B. Downregulated Genes',
                         fontsize=14, fontweight='bold', pad=20, loc='left')
        ax2.grid(True, alpha=0.3, linewidth=0.8, axis='x')
        ax2.set_axisbelow(True)

        # Add fold change values on bars
        for i, (bar, fc, pval) in enumerate(zip(bars2, top_down['log2FC'], top_down['padj'])):
            width = bar.get_width()
            ax2.text(width - 0.05, bar.get_y() + bar.get_height() / 2,
                         f'{fc:.2f}', va='center', ha='right', fontsize=9, fontweight='bold')

        # Set x-axis limits with padding
        ax2.set_xlim(min(top_down['log2FC']) * 1.3, 0)

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

            # Use axis transform to precisely position stars
            transform = ax2.get_yaxis_transform()
            ax2.text(1.02, i, sig_text, transform=transform, va='center', ha='left',
                         fontsize=10, fontweight='bold', color=colors_pub['blue'])

    # Add significance legend below the plots
    from matplotlib.lines import Line2D
    sig_legend_elements = [
        Line2D([0], [0], color='w', label='***  p < 0.001'),
        Line2D([0], [0], color='w', label='**  p < 0.01'),
        Line2D([0], [0], color='w', label='*  p < 0.05'),
    ]

    fig.legend(handles=sig_legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.01),
               ncol=3, frameon=False, title="Significance levels:")

    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.subplots_adjust(left=0.2, wspace=0.8)

    output_path = output_dir / 'top_genes_barplot'
    print(f"  Saving top genes bar plot to {output_path}...")
    plt.savefig(f'{output_path}.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.eps', bbox_inches='tight', facecolor='white')
    # plt.show()
    plt.close()


def create_summary_stats_plot(deg_df: pd.DataFrame, significant_genes: pd.DataFrame, significant_corr: pd.DataFrame,
                              output_dir: Path):
    """
    Create analysis results summary plot, styled for publication.
    """
    print("Creating analysis results summary plot...")

    # --- Data Preparation ---
    up_genes = len(significant_genes[significant_genes['log2FC'] > 0])
    down_genes = len(significant_genes[significant_genes['log2FC'] < 0])
    total_degs = up_genes + down_genes

    pos_corr = len(significant_corr[significant_corr['rho'] > 0])
    neg_corr = len(significant_corr[significant_corr['rho'] < 0])
    total_corr = pos_corr + neg_corr

    # --- Plotting ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Summary of Differential Expression and Correlation Analysis', fontsize=18, fontweight='bold')

    # --- 1. DEG Statistics (Replaced Pie with Bar Chart) ---
    ax = axes[0, 0]
    deg_labels = ['Upregulated', 'Downregulated']
    deg_counts = [up_genes, down_genes]
    deg_colors = [colors_pub['red'], colors_pub['blue']]
    bars = ax.barh(deg_labels, deg_counts, color=deg_colors, height=0.6, edgecolor='black', linewidth=1.2)
    ax.set_title('A. Significant DEGs Distribution', fontsize=14, fontweight='bold', loc='left')
    ax.set_xlabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.tick_params(axis='y', labelsize=12)

    # Add annotations
    if total_degs > 0:
        for i, bar in enumerate(bars):
            width = bar.get_width()
            percentage = width / total_degs * 100
            ax.text(width + max(deg_counts)*0.01, bar.get_y() + bar.get_height()/2, f'{width} ({percentage:.1f}%)',
                    va='center', ha='left', fontsize=10, fontweight='bold')
    
    ax.set_xlim(0, max(deg_counts) * 1.25)
    sns.despine(ax=ax)
    ax.grid(axis='x', linestyle='--', alpha=0.6)

    # --- 2. Effect size distribution ---
    ax = axes[0, 1]
    fc_bins = [0, 0.5, 1, 1.5, 2, float('inf')]
    fc_labels = ['0-0.5', '0.5-1', '1-1.5', '1.5-2', '>2']
    fc_counts = []
    if not significant_genes.empty:
    for i in range(len(fc_bins) - 1):
        count = len(significant_genes[(abs(significant_genes['log2FC']) >= fc_bins[i]) &
                                      (abs(significant_genes['log2FC']) < fc_bins[i + 1])])
        fc_counts.append(count)
    else:
        fc_counts = [0] * len(fc_labels)

    bars = ax.bar(fc_labels, fc_counts, color=colors_pub['green'], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.set_title('B. Effect Size Distribution of DEGs', fontsize=14, fontweight='bold', loc='left')
    ax.set_xlabel('|log2(Fold Change)|', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.tick_params(axis='x', rotation=0, labelsize=12)

    # Add annotations on bars
    if sum(fc_counts) > 0:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax.set_ylim(0, max(fc_counts) * 1.15)
    sns.despine(ax=ax)
    ax.grid(axis='y', linestyle='--', alpha=0.6)

    # --- 3. Significance level distribution of DEGs ---
    ax = axes[1, 0]
    p_bins = [0, 0.001, 0.01, 0.05]
    p_labels = ['< 0.001', '0.001-0.01', '0.01-0.05']
    p_counts = []
    if not significant_genes.empty:
    for i in range(len(p_bins) - 1):
            count = len(significant_genes[(significant_genes['padj'] >= p_bins[i]) &
                                          (significant_genes['padj'] < p_bins[i + 1])])
        p_counts.append(count)
    else:
        p_counts = [0] * len(p_labels)
        
    bars = ax.bar(p_labels, p_counts, color=colors_pub['gray'], alpha=0.8, edgecolor='black', linewidth=1.2)
    ax.set_title('C. Significance Level of DEGs', fontsize=14, fontweight='bold', loc='left')
    ax.set_xlabel('Adjusted P-value', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.tick_params(axis='x', rotation=0, labelsize=12)

    if sum(p_counts) > 0:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height}',
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax.set_ylim(0, max(p_counts) * 1.15)
    sns.despine(ax=ax)
    ax.grid(axis='y', linestyle='--', alpha=0.6)

    # --- 4. Correlation statistics (Replaced Pie with Bar Chart) ---
    ax = axes[1, 1]
    corr_labels = ['Positive', 'Negative']
    corr_counts = [pos_corr, neg_corr]
    corr_colors = [colors_pub['orange'], colors_pub['purple']]
    bars = ax.barh(corr_labels, corr_counts, color=corr_colors, height=0.6, edgecolor='black', linewidth=1.2)
    ax.set_title('D. Significant Correlations Distribution', fontsize=14, fontweight='bold', loc='left')
    ax.set_xlabel('Number of Genes', fontsize=12, fontweight='bold')
    ax.tick_params(axis='y', labelsize=12)

    if total_corr > 0:
        for i, bar in enumerate(bars):
            width = bar.get_width()
            percentage = width / total_corr * 100
            ax.text(width + max(corr_counts)*0.01, bar.get_y() + bar.get_height()/2, f'{width} ({percentage:.1f}%)',
                    va='center', ha='left', fontsize=10, fontweight='bold')
    
    ax.set_xlim(0, max(corr_counts) * 1.25)
    sns.despine(ax=ax)
    ax.grid(axis='x', linestyle='--', alpha=0.6)

    # --- Final Layout Adjustment ---
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    output_path = output_dir / 'summary_statistics'
    print(f"  Saving summary statistics plot to {output_path}...")
    plt.savefig(f'{output_path}.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight', facecolor='white')
    plt.savefig(f'{output_path}.eps', bbox_inches='tight', facecolor='white')
    # plt.show()
    plt.close()


def main():
    """
    Main visualization workflow
    """
    parser = argparse.ArgumentParser(description="Visualize results from BRM gene RNA-seq analysis.")
    parser.add_argument(
        "--input_dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "brm_analysis_csv",
        help="Directory containing the analysis result CSV files. Default: results/brm_analysis_csv"
    )
    parser.add_argument(
        "--expression_file",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "working.csv",
        help="Path to the processed expression data file. Default: data/processed/working.csv"
    )
    parser.add_argument(
        "--target_gene",
        type=str,
        default="SMARCA2",
        help="The target gene symbol (e.g., SMARCA2)."
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "visualize_results",
        help="Directory to save the generated plots. Default: results/visualize_results"
    )
    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=== BRM Gene RNA-seq Data Analysis Results Visualization ===\n")
    print(f"Configuration:")
    print(f"  - Input Directory : {args.input_dir.relative_to(PROJECT_ROOT)}")
    print(f"  - Output Directory: {args.output_dir.relative_to(PROJECT_ROOT)}")
    print(f"  - Expression Data : {args.expression_file.relative_to(PROJECT_ROOT)}")
    print(f"  - Target Gene     : {args.target_gene}\n")

    # Load data
    (
        deg_df,
        significant_genes,
        corr_df,
        significant_corr,
        expr_filtered,
        target_gene_full_name
    ) = load_data(
        input_dir=args.input_dir,
        expression_file=args.expression_file,
        target_gene=args.target_gene
    )

    if expr_filtered.empty or target_gene_full_name is None:
        print("\nCould not proceed with visualization due to missing data or target gene.")
        return

    # Generate various plots
    print("\n1. Creating volcano plot...")
    create_volcano_plot(deg_df, args.output_dir, args.target_gene)

    print("\n2. Creating expression distribution plot...")
    create_expression_distribution(expr_filtered, target_gene_full_name, args.output_dir, args.target_gene)

    print("\n3. Creating DEGs heatmap...")
    create_deg_heatmap(expr_filtered, significant_genes, target_gene_full_name, args.output_dir, args.target_gene,
                       top_n=50)

    print("\n4. Creating correlation analysis plots...")
    create_correlation_plots(corr_df, significant_corr, args.output_dir, args.target_gene)

    print("\n5. Creating top genes bar plot...")
    create_top_genes_barplot(significant_genes, args.output_dir, top_n=15)

    print("\n6. Creating summary statistics plot...")
    create_summary_stats_plot(deg_df, significant_genes, significant_corr, args.output_dir)

    print("\n=== Visualization Complete ===")
    print(f"All plots have been saved to '{args.output_dir.relative_to(PROJECT_ROOT)}' directory.")
    print("Generated plots include:")
    print("- volcano_plot")
    print("- smarca2_expression_distribution")
    print("- deg_heatmap")
    print("- correlation_analysis")
    print("- top_genes_barplot")
    print("- summary_statistics")


if __name__ == "__main__":
    main() 