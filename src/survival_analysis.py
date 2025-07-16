
import argparse
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

# 忽略代码运行过程中可能出现的警告信息
warnings.filterwarnings('ignore')

# 定义项目的根目录，即当前脚本所在位置的上两级目录
PROJECT_ROOT = Path(__file__).resolve().parents[1]

def perform_survival_analysis(expression_file: Path, clinical_file: Path, target_gene: str, output_dir: Path):
    """
    根据基因表达水平执行生存分析。

    Args:
        expression_file (Path): 基因表达数据文件的路径 (CSV格式)。
        clinical_file (Path): 临床数据文件的路径。
        target_gene (str): 要分析的目标基因符号 (例如, 'SMARCA2')。
        output_dir (Path): 保存图表的输出目录。
    """
    print("开始生存分析...")

    # 1. 加载数据
    print("  - 正在加载基因表达和临床数据...")
    try:
        # 加载基因表达矩阵，索引列为基因名
        expr_df = pd.read_csv(expression_file, index_col=0)
        # 加载临床数据，注意它是制表符分隔的，并且需要转置（行是样本，列是临床特征）
        clin_df = pd.read_csv(clinical_file, sep='\t', index_col=0).T
        print(f"  - 原始表达数据维度: {expr_df.shape}")
        print(f"  - 原始临床数据维度: {clin_df.shape}")
    except FileNotFoundError as e:
        print(f"数据加载错误: {e}")
        return

    # 2. 在表达数据中查找目标基因的完整名称
    target_gene_full_name = None
    for gene in expr_df.index:
        # 匹配基因符号的开头部分，忽略大小写
        if str(gene).upper().startswith(target_gene.upper()):
            target_gene_full_name = gene
            print(f"  - 在表达数据中找到目标基因: {target_gene_full_name}")
            break
    
    if not target_gene_full_name:
        print(f"Error: Target gene '{target_gene}' not found in expression data.")
        return
        
    # 3. 准备数据以便合并 (Standardize Sample IDs)
    print("  - Standardizing sample IDs for merging...")
    # 提取目标基因的表达数据
    target_expr = expr_df.loc[target_gene_full_name]
    target_expr.name = f"{target_gene}_expression"
    
    # 标准化表达数据的样本ID：大写，替换点，截取前12位
    target_expr.index = target_expr.index.str.upper().str.replace('.', '-').str.slice(0, 12)
    
    # 标准化临床数据的样本ID
    clin_df.index = clin_df.index.str.upper().str.slice(0, 12)
    print("  - 示例标准化ID (表达):", target_expr.index[:3].tolist())
    print("  - 示例标准化ID (临床):", clin_df.index[:3].tolist())

    # 在合并前移除因ID截断可能产生的重复项
    original_expr_count = len(target_expr)
    original_clin_count = len(clin_df)
    target_expr = target_expr[~target_expr.index.duplicated(keep='first')]
    clin_df = clin_df[~clin_df.index.duplicated(keep='first')]
    print(f"  - 表达数据去重后: {len(target_expr)} / {original_expr_count} 个唯一样本")
    print(f"  - 临床数据去重后: {len(clin_df)} / {original_clin_count} 个唯一样本")
    
    # 4. Merge DataFrames
    print("  - 正在合并表达和临床数据...")
    # 基于样本ID（索引）将两个DataFrame合并
    merged_df = pd.merge(target_expr, clin_df, left_index=True, right_index=True)
    print(f"  - 合并后数据维度: {merged_df.shape}")
    
    # 在合并后立即进行检查
    print(f"  - Found {len(merged_df)} common samples after merging.")
    if merged_df.empty:
        print("错误: 在表达数据和临床数据之间没有找到共同的样本。请检查文件和ID格式。")
        return

    # 5. 准备生存分析所需的列
    print("  - 正在处理生存数据...")
    # 'event'列：事件状态。'vital_status'中1代表死亡，0代表存活。
    merged_df['event'] = pd.to_numeric(merged_df['vital_status'], errors='coerce')
    
    # 'time'列：生存时间。需要结合'days_to_death'和'days_to_last_followup'
    days_to_death = pd.to_numeric(merged_df['days_to_death'], errors='coerce')
    days_to_last_followup = pd.to_numeric(merged_df['days_to_last_followup'], errors='coerce')
    
    # 如果患者已死亡(event=1)，生存时间是'days_to_death'
    # 如果患者存活(event=0)，生存时间是'days_to_last_followup'
    merged_df['time'] = np.where(merged_df['event'] == 1, days_to_death, days_to_last_followup)
    
    # 记录处理前的样本数
    samples_before_dropna = len(merged_df)
    # 丢弃没有有效时间和事件信息的行
    merged_df.dropna(subset=['time', 'event'], inplace=True)
    samples_after_dropna = len(merged_df)
    print(f"  - 因生存数据缺失，丢弃了 {samples_before_dropna - samples_after_dropna} 个样本。")
    
    if merged_df.empty:
        print("错误: 清理后没有可用的生存数据。")
        return
        
    print(f"  - 找到 {len(merged_df)} 个具有完整生存数据的样本。")

    # 6. 根据基因表达的四分位数创建高表达组和低表达组
    low_quantile = merged_df[f"{target_gene}_expression"].quantile(0.25)
    high_quantile = merged_df[f"{target_gene}_expression"].quantile(0.75)
    
    # 根据分位数定义分组
    merged_df['group'] = np.nan
    merged_df.loc[merged_df[f"{target_gene}_expression"] >= high_quantile, 'group'] = 'High'
    merged_df.loc[merged_df[f"{target_gene}_expression"] <= low_quantile, 'group'] = 'Low'
    
    print(f"  - 已根据表达值的四分位数进行分组 (低-25%: <= {low_quantile:.2f}, 高-25%: >= {high_quantile:.2f})")

    # 筛选出高表达和低表达组的样本用于分析
    analysis_data = merged_df.dropna(subset=['group'])
    
    print(f"  - 原始样本数 (有生存数据): {len(merged_df)}")
    print(f"  - 分析用样本数 (高/低表达组): {len(analysis_data)}")
    high_group_count = len(analysis_data[analysis_data['group'] == 'High'])
    low_group_count = len(analysis_data[analysis_data['group'] == 'Low'])
    print(f"  - 高表达组 (top 25%) 样本数: {high_group_count}")
    print(f"  - 低表达组 (bottom 25%) 样本数: {low_group_count}")

    if high_group_count == 0 or low_group_count == 0:
        print("错误: 高表达组或低表达组中没有样本，无法进行分析。")
        return
    
    # 7. 执行生存分析
    print("  - 正在拟合Kaplan-Meier生存曲线并执行log-rank检验...")
    
    # 创建KaplanMeierFitter实例
    kmf = KaplanMeierFitter()
    
    # 分别准备高表达组和低表达组的数据
    high_group = analysis_data[analysis_data['group'] == 'High']
    low_group = analysis_data[analysis_data['group'] == 'Low']
    
    # 拟合高表达组数据并绘制
    kmf.fit(durations=high_group['time'], event_observed=high_group['event'], label=f'High {target_gene} (n={high_group_count})')
    ax = kmf.plot(figsize=(10, 8), ci_show=False) # ci_show=False 不显示置信区间
    
    # 在同一张图上拟合低表达组数据并绘制
    kmf.fit(durations=low_group['time'], event_observed=low_group['event'], label=f'Low {target_gene} (n={low_group_count})')
    kmf.plot(ax=ax, ci_show=False)
    
    # 执行log-rank检验，比较两组生存曲线的差异
    results = logrank_test(
        durations_A=high_group['time'], event_observed_A=high_group['event'],
        durations_B=low_group['time'], event_observed_B=low_group['event']
    )
    
    # 8. 美化并保存图表
    print("  - 正在生成图表...")
    plt.title(f"Survival Analysis based on {target_gene} Expression in OV", fontsize=16, fontweight='bold')
    plt.xlabel("Time (days)", fontsize=12, fontweight='bold')
    plt.ylabel("Survival Probability", fontsize=12, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # 在图表上标注log-rank检验的p值
    plt.text(0.5, 0.1, f'Log-rank p-value: {results.p_value:.4f}', 
             transform=ax.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.8))
             
    plt.legend(fancybox=True, shadow=True)
    
    # 定义输出路径并保存多种格式的图片
    output_path = output_dir / f"{target_gene}_survival_analysis"
    print(f"  - 正在保存生存分析图表至 {output_path.name}.png")
    plt.savefig(f"{output_path}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}.pdf", bbox_inches='tight')
    plt.show()
    plt.close()

    # 打印log-rank检验的详细结果
    print(f"\nLog-rank检验结果:\n{results.summary}")


def main():
    """
    生存分析主流程。
    """
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="Perform survival analysis based on gene expression.")
    parser.add_argument(
        "--expression_file",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "working.csv",
        help="处理后的基因表达数据文件路径。"
    )
    parser.add_argument(
        "--clinical_file",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "OV.clin.merged.picked.txt",
        help="临床数据文件路径。"
    )
    parser.add_argument(
        "--target_gene",
        type=str,
        default="SMARCA2",
        help="目标基因的符号 (例如, SMARCA2)。"
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "survival_analysis",
        help="保存分析图表的目录。"
    )
    args = parser.parse_args()

    # 创建输出目录，如果不存在的话
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=== 生存分析工作流 ===\n")
    # 调用核心函数执行分析
    perform_survival_analysis(
        expression_file=args.expression_file,
        clinical_file=args.clinical_file,
        target_gene=args.target_gene,
        output_dir=args.output_dir
    )
    print("\n=== 生存分析完成 ===")


if __name__ == "__main__":
    # 当该脚本被直接执行时，运行main函数
    main() 