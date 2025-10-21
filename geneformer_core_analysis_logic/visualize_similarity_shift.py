#!/usr/bin/env python3
"""
可视化相似度变化分析结果
分析基因敲除是否能将耐药细胞状态推向敏感细胞状态
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

def setup_logging(results_dir):
    """设置日志"""
    log_file = results_dir / f"visualization_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    def log_func(message):
        timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
        log_message = f"{timestamp} {message}"
        print(log_message)
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(log_message + '\n')
    
    return log_func

def load_and_analyze_data(csv_path, log_func):
    """加载和分析数据"""
    log_func(f"加载数据: {csv_path}")
    df = pd.read_csv(csv_path)
    
    log_func(f"数据形状: {df.shape}")
    log_func(f"包含基因: {df['gene_name'].unique()}")
    
    # 计算每个基因的统计信息
    gene_stats = []
    for gene in df['gene_name'].unique():
        gene_data = df[df['gene_name'] == gene]
        
        stats = {
            'gene_name': gene,
            'gene_id': gene_data['gene_id'].iloc[0],
            'cell_count': len(gene_data),
            'mean_shift': gene_data['similarity_shift'].mean(),
            'median_shift': gene_data['similarity_shift'].median(),
            'std_shift': gene_data['similarity_shift'].std(),
            'positive_shift_count': (gene_data['similarity_shift'] > 0).sum(),
            'negative_shift_count': (gene_data['similarity_shift'] < 0).sum(),
            'positive_ratio': (gene_data['similarity_shift'] > 0).mean(),
            'min_shift': gene_data['similarity_shift'].min(),
            'max_shift': gene_data['similarity_shift'].max()
        }
        gene_stats.append(stats)
    
    gene_stats_df = pd.DataFrame(gene_stats)
    log_func("\n基因统计信息:")
    for _, row in gene_stats_df.iterrows():
        log_func(f"{row['gene_name']}: 细胞数={row['cell_count']}, "
                f"平均变化={row['mean_shift']:.6f}, "
                f"正向变化比例={row['positive_ratio']:.2%}")
    
    return df, gene_stats_df

def create_boxplot_visualization(df, output_dir, log_func):
    """创建箱线图可视化"""
    log_func("Creating boxplot visualization...")
    
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.unicode_minus'] = False
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
    
    # Boxplot 1: Similarity Shift
    box_data1 = [df[df['gene_name'] == gene]['similarity_shift'].values 
                 for gene in df['gene_name'].unique()]
    
    bp1 = ax1.boxplot(box_data1, labels=df['gene_name'].unique(), patch_artist=True,
                      showfliers=False) # Hide outliers for clarity
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(bp1['boxes'])))
    for patch, color in zip(bp1['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        
    for median in bp1['medians']:
        median.set_color('black')
    
    ax1.set_xlabel('Gene', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Similarity Shift (Post-KO - Pre-KO)', fontsize=14, fontweight='bold')
    ax1.set_title('Effect of Gene Knockout on Cell Similarity', fontsize=16, fontweight='bold')
    ax1.axhline(y=0, color='red', linestyle='--', alpha=0.9, label='No Change')
    ax1.legend()
    
    plt.setp(ax1.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    plt.setp(ax1.get_yticklabels(), fontsize=12)
    
    # Boxplot 2: Post-Knockout Absolute Similarity
    box_data2 = [df[df['gene_name'] == gene]['post_knockout_similarity'].values 
                 for gene in df['gene_name'].unique()]
    
    bp2 = ax2.boxplot(box_data2, labels=df['gene_name'].unique(), patch_artist=True,
                      showfliers=False)
    
    for patch, color in zip(bp2['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)

    for median in bp2['medians']:
        median.set_color('black')
        
    ax2.set_xlabel('Gene', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Post-KO Similarity to Sensitive Group', fontsize=14, fontweight='bold')
    ax2.set_title('Post-Knockout Similarity Distribution', fontsize=16, fontweight='bold')
    
    plt.setp(ax2.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    plt.setp(ax2.get_yticklabels(), fontsize=12)
    
    plt.tight_layout(pad=3.0)
    
    output_path = output_dir / 'similarity_shift_boxplots.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    log_func(f"Boxplot saved to: {output_path}")

def create_violin_plot(df, output_dir, log_func):
    """创建小提琴图"""
    log_func("Creating violin plot...")
    
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['font.family'] = 'Arial'
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # 创建小提琴图
    sns.violinplot(data=df, x='gene_name', y='similarity_shift', ax=ax,
                   palette='viridis', inner='quartile', linewidth=1.5)
    
    ax.set_xlabel('Gene', fontsize=14, fontweight='bold')
    ax.set_ylabel('Similarity Shift', fontsize=14, fontweight='bold')
    ax.set_title('Distribution of Similarity Shift upon Gene Knockout', fontsize=16, fontweight='bold')
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.9)
    
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    plt.setp(ax.get_yticklabels(), fontsize=12)
    
    plt.tight_layout()
    
    output_path = output_dir / 'similarity_shift_violinplot.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    log_func(f"Violin plot saved to: {output_path}")

def create_heatmap(df, output_dir, log_func):
    """创建热力图显示基因效应"""
    log_func("Creating gene effect heatmap...")
    
    # 定义效应类别和阈值
    effect_categories = {
        'Strong Positive (>0.001)': (0.001, np.inf),
        'Moderate Positive (0.0001-0.001)': (0.0001, 0.001),
        'Weak Positive (0-0.0001)': (0, 0.0001),
        'Weak Negative (-0.0001-0)': (-0.0001, 0),
        'Moderate Negative (-0.001--0.0001)': (-0.001, -0.0001),
        'Strong Negative (<-0.001)': (-np.inf, -0.001)
    }

    gene_effects = []
    for gene in df['gene_name'].unique():
        gene_data = df[df['gene_name'] == gene]
        effects = {'Gene': gene, 'Total Cells': len(gene_data)}
        for name, (lower, upper) in effect_categories.items():
            effects[name] = ((gene_data['similarity_shift'] > lower) & 
                             (gene_data['similarity_shift'] <= upper)).sum()
        gene_effects.append(effects)
    
    effects_df = pd.DataFrame(gene_effects)
    
    # 转换为比例
    for col in effect_categories.keys():
        effects_df[f'{col} Ratio'] = effects_df[col] / effects_df['Total Cells']
    
    heatmap_data = effects_df[[f'{cat} Ratio' for cat in effect_categories.keys()]].T
    heatmap_data.columns = effects_df['Gene']
    heatmap_data.index = list(effect_categories.keys())
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    sns.heatmap(heatmap_data, annot=True, fmt='.1%', cmap='coolwarm', 
                center=0, ax=ax, cbar_kws={'label': 'Percentage of Cells'},
                linewidths=.5)
    
    ax.set_title('Heatmap of Gene Knockout Effect Distribution', fontsize=16, fontweight='bold')
    ax.set_xlabel('Gene', fontsize=14, fontweight='bold')
    ax.set_ylabel('Effect Category', fontsize=14, fontweight='bold')
    
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=12)
    
    plt.tight_layout()
    
    output_path = output_dir / 'gene_effect_heatmap.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    log_func(f"Heatmap saved to: {output_path}")
    
    effects_df.to_csv(output_dir / 'gene_effect_statistics.csv', index=False)
    log_func(f"Effect statistics saved to: {output_dir / 'gene_effect_statistics.csv'}")

def create_scatter_plots(df, output_dir, log_func):
    """创建散点图"""
    log_func("Creating scatter plots...")
    
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['font.family'] = 'Arial'
    
    fig = plt.figure(figsize=(16, 14))
    gs = fig.add_gridspec(2, 2)

    # Scatter Plot 1: Pre-KO vs Post-KO Similarity
    ax1 = fig.add_subplot(gs[0, 0])
    sns.scatterplot(data=df, x='pre_knockout_similarity', y='post_knockout_similarity',
                    hue='gene_name', palette='viridis', alpha=0.7, s=50, ax=ax1)
    ax1.plot([df['pre_knockout_similarity'].min(), 1], 
             [df['pre_knockout_similarity'].min(), 1], 
             'r--', alpha=0.8, label='y=x (No Change)')
    ax1.set_xlabel('Pre-KO Similarity', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Post-KO Similarity', fontsize=14, fontweight='bold')
    ax1.set_title('Pre- vs Post-Knockout Similarity', fontsize=16, fontweight='bold')
    ax1.legend(title='Gene', fontsize=10)
    
    # Scatter Plot 2: Baseline Similarity vs. Shift
    ax2 = fig.add_subplot(gs[0, 1])
    sns.scatterplot(data=df, x='pre_knockout_similarity', y='similarity_shift',
                    hue='gene_name', palette='viridis', alpha=0.7, s=50, ax=ax2)
    ax2.axhline(y=0, color='red', linestyle='--', alpha=0.8)
    ax2.set_xlabel('Pre-KO Similarity', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Similarity Shift', fontsize=14, fontweight='bold')
    ax2.set_title('Baseline Similarity vs. Effect Size', fontsize=16, fontweight='bold')
    ax2.legend(title='Gene', fontsize=10)
    
    # Scatter Plot 3: Mean Effect Size per Gene
    ax3 = fig.add_subplot(gs[1, 0])
    gene_means = df.groupby('gene_name')['similarity_shift'].agg(['mean', 'std']).reset_index()
    
    ax3.errorbar(gene_means['gene_name'], gene_means['mean'], yerr=gene_means['std'], 
                 fmt='o', capsize=5, capthick=2, markersize=8, color='darkblue')
    ax3.axhline(y=0, color='red', linestyle='--', alpha=0.8)
    ax3.set_ylabel('Mean Similarity Shift ± SD', fontsize=14, fontweight='bold')
    ax3.set_title('Mean Knockout Effect per Gene', fontsize=16, fontweight='bold')
    plt.setp(ax3.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    
    # Scatter Plot 4: Cell Count vs. Mean Effect
    ax4 = fig.add_subplot(gs[1, 1])
    gene_summary = df.groupby('gene_name').agg({
        'similarity_shift': ['count', 'mean']
    }).reset_index()
    gene_summary.columns = ['gene_name', 'cell_count', 'mean_shift']
    
    sns.scatterplot(data=gene_summary, x='cell_count', y='mean_shift', 
                    size='cell_count', hue='gene_name', sizes=(100, 500), 
                    palette='viridis', ax=ax4, legend=False)

    for i, row in gene_summary.iterrows():
        ax4.annotate(row['gene_name'],
                     xy=(row['cell_count'], row['mean_shift']),
                     xytext=(5, 0),
                     textcoords='offset points',
                     ha='left',
                     va='center',
                     fontsize=10)
    
    ax4.axhline(y=0, color='red', linestyle='--', alpha=0.8)
    ax4.set_xlabel('Number of Cells Expressing Gene', fontsize=14, fontweight='bold')
    ax4.set_ylabel('Mean Similarity Shift', fontsize=14, fontweight='bold')
    ax4.set_title('Cell Count vs. Mean Effect', fontsize=16, fontweight='bold')

    plt.tight_layout(pad=3.0)
    
    output_path = output_dir / 'similarity_shift_scatterplots.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    log_func(f"Scatter plots saved to: {output_path}")

def generate_summary_report(df, gene_stats_df, output_dir, log_func):
    """生成总结报告"""
    log_func("Generating summary report...")
    
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("Analysis Report: Effect of Gene Knockout on Cell State Transition")
    report_lines.append("=" * 80)
    report_lines.append(f"Analysis Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"Total Data Points Analyzed: {len(df)}")
    report_lines.append(f"Number of Genes Analyzed: {len(df['gene_name'].unique())}")
    report_lines.append("")
    
    report_lines.append("1. Overview of Gene Effects:")
    report_lines.append("-" * 50)
    for _, row in gene_stats_df.iterrows():
        effect_direction = "Promotes Transition" if row['mean_shift'] > 0 else "Inhibits Transition"
        report_lines.append(f"Gene: {row['gene_name']} (ID: {row['gene_id']})")
        report_lines.append(f"  - Cells Analyzed: {row['cell_count']}")
        report_lines.append(f"  - Mean Effect (Similarity Shift): {row['mean_shift']:.6f} ({effect_direction})")
        report_lines.append(f"  - Proportion of Cells with Positive Shift: {row['positive_ratio']:.2%}")
        report_lines.append(f"  - Effect Range: [{row['min_shift']:.6f}, {row['max_shift']:.6f}]")
        report_lines.append("")
    
    report_lines.append("2. Key Findings:")
    report_lines.append("-" * 50)
    
    best_gene = gene_stats_df.loc[gene_stats_df['mean_shift'].idxmax()]
    worst_gene = gene_stats_df.loc[gene_stats_df['mean_shift'].idxmin()]
    
    report_lines.append(f"Gene with Strongest Positive Effect: {best_gene['gene_name']} (Mean Shift: {best_gene['mean_shift']:.6f})")
    report_lines.append(f"Gene with Strongest Negative Effect: {worst_gene['gene_name']} (Mean Shift: {worst_gene['mean_shift']:.6f})")
    
    strong_positive = (df['similarity_shift'] > 0.001).sum()
    strong_negative = (df['similarity_shift'] < -0.001).sum()
    
    report_lines.append(f"Cells with Strong Positive Shift (>0.001): {strong_positive} ({strong_positive/len(df):.2%})")
    report_lines.append(f"Cells with Strong Negative Shift (<-0.001): {strong_negative} ({strong_negative/len(df):.2%})")
    report_lines.append("")
    
    report_lines.append("3. Statistical Summary:")
    report_lines.append("-" * 50)
    report_lines.append(f"Overall Mean Similarity Shift: {df['similarity_shift'].mean():.6f}")
    report_lines.append(f"Overall Median Similarity Shift: {df['similarity_shift'].median():.6f}")
    report_lines.append(f"Standard Deviation of Shift: {df['similarity_shift'].std():.6f}")
    report_lines.append(f"Overall Proportion of Positive Shifts: {(df['similarity_shift'] > 0).mean():.2%}")
    report_lines.append("")
    
    report_lines.append("4. Conclusion:")
    report_lines.append("-" * 50)
    if df['similarity_shift'].mean() > 0:
        report_lines.append("On average, knockout of the targeted genes tends to increase the similarity of refractory cells to the sensitive state.")
        report_lines.append("This suggests that inhibiting these genes could be a potential strategy to overcome drug resistance.")
    else:
        report_lines.append("On average, knockout of the targeted genes tends to decrease the similarity of refractory cells to the sensitive state.")
        report_lines.append("This suggests that inhibiting these genes may not be a viable strategy for overcoming drug resistance in this context.")
    
    report_lines.append("")
    report_lines.append("=" * 80)
    
    report_path = output_dir / 'analysis_summary_report.txt'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    log_func(f"Summary report saved to: {report_path}")
    
    # 打印关键发现到控制台
    log_func("\n" + "="*20 + " KEY FINDINGS " + "="*20)
    log_func(f"Gene with Strongest Positive Effect: {best_gene['gene_name']} (Mean Shift: {best_gene['mean_shift']:.6f})")
    log_func(f"Gene with Strongest Negative Effect: {worst_gene['gene_name']} (Mean Shift: {worst_gene['mean_shift']:.6f})")
    log_func(f"Overall Mean Similarity Shift: {df['similarity_shift'].mean():.6f}")
    log_func("=" * 52)

def main():
    # Setup paths
    project_root = Path(__file__).resolve().parent.parent
    results_dir = project_root / 'results' / 'refractory_to_sensitive_shift_analysis'
    
    viz_dir = results_dir / 'visualizations'
    viz_dir.mkdir(exist_ok=True)
    
    log_func = setup_logging(viz_dir)
    
    log_func("=" * 60)
    log_func("Starting Visualization of Similarity Shift Analysis")
    log_func("=" * 60)
    
    csv_path = results_dir / 'similarity_shift_results.csv'
    df, gene_stats_df = load_and_analyze_data(csv_path, log_func)
    
    # Create visualizations
    create_boxplot_visualization(df, viz_dir, log_func)
    create_violin_plot(df, viz_dir, log_func)
    create_heatmap(df, viz_dir, log_func)
    create_scatter_plots(df, viz_dir, log_func)
    
    generate_summary_report(df, gene_stats_df, viz_dir, log_func)
    
    log_func("=" * 60)
    log_func("Visualization analysis complete!")
    log_func(f"All results saved to: {viz_dir}")
    log_func("=" * 60)

if __name__ == "__main__":
    main()
