
import argparse
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt

# 忽略代码运行过程中可能出现的警告信息
warnings.filterwarnings('ignore')

# 定义项目的根目录
PROJECT_ROOT = Path(__file__).resolve().parents[1]

def perform_cox_regression(expression_file: Path, clinical_file: Path, target_gene: str, output_dir: Path):
    """
    根据基因表达水平和临床变量执行Cox比例风险回归分析。

    Args:
        expression_file (Path): 基因表达数据文件的路径 (CSV格式)。
        clinical_file (Path): 临床数据文件的路径。
        target_gene (str): 要分析的目标基因符号 (例如, 'SMARCA2')。
        output_dir (Path): 保存图表的输出目录。
    """
    print("=== 开始Cox比例风险回归分析 ===")

    # 1. 加载数据
    print("\n[1. 数据加载]")
    try:
        expr_df = pd.read_csv(expression_file, index_col=0)
        clin_df = pd.read_csv(clinical_file, sep='\t', index_col=0).T
        print(f"  - 原始表达数据维度: {expr_df.shape}")
        print(f"  - 原始临床数据维度: {clin_df.shape}")
    except FileNotFoundError as e:
        print(f"  - 数据加载错误: {e}")
        return

    # 2. 准备和合并数据
    print("\n[2. 数据准备与合并]")
    # 查找目标基因
    target_gene_full_name = next((gene for gene in expr_df.index if str(gene).upper().startswith(target_gene.upper())), None)
    if not target_gene_full_name:
        print(f"  - 错误: 在表达数据中未找到目标基因 '{target_gene}'。")
        return
    print(f"  - 在表达数据中找到目标基因: {target_gene_full_name}")

    # 标准化样本ID并合并
    target_expr = expr_df.loc[target_gene_full_name]
    target_expr.name = f"{target_gene}_expression"
    target_expr.index = target_expr.index.str.upper().str.replace('.', '-').str.slice(0, 12)
    clin_df.index = clin_df.index.str.upper().str.slice(0, 12)
    
    target_expr = target_expr[~target_expr.index.duplicated(keep='first')]
    clin_df = clin_df[~clin_df.index.duplicated(keep='first')]
    
    merged_df = pd.merge(pd.DataFrame(target_expr), clin_df, left_index=True, right_index=True)
    print(f"  - 找到 {len(merged_df)} 个共同样本进行合并。")

    if merged_df.empty:
        print("  - 错误: 表达和临床数据之间无共同样本。")
        return

    # 3. 特征工程与数据清洗
    print("\n[3. 特征工程与数据清洗]")
    # 创建生存时间和事件列
    merged_df['event'] = pd.to_numeric(merged_df['vital_status'], errors='coerce')
    days_to_death = pd.to_numeric(merged_df['days_to_death'], errors='coerce')
    days_to_last_followup = pd.to_numeric(merged_df['days_to_last_followup'], errors='coerce')
    merged_df['time'] = np.where(merged_df['event'] == 1, days_to_death, days_to_last_followup)

    # 计算年龄
    current_year = pd.to_datetime('today').year
    merged_df['age_at_diagnosis'] = pd.to_numeric(merged_df['date_of_initial_pathologic_diagnosis'], errors='coerce') - pd.to_numeric(merged_df['years_to_birth'], errors='coerce')
    print(f"  - 已计算'age_at_diagnosis' (诊断时年龄)。")
    
    # 由于 'pathologic_stage' 数据缺失严重，本次分析仅使用年龄作为协变量
    print(f"  - 注意: 'pathologic_stage' 数据缺失，将仅使用年龄进行多变量校正。")

    print("\n  - 正在检查关键列中的缺失值...")
    key_cols_for_check = ['time', 'event', 'age_at_diagnosis', f"{target_gene}_expression"]
    for col in key_cols_for_check:
        if col in merged_df.columns:
            print(f"    - '{col}' 缺失值: {merged_df[col].isnull().sum()} / {len(merged_df)}")
        else:
            print(f"    - '{col}' 列不存在。")

    # 确定用于模型的列
    cox_columns = [
        f"{target_gene}_expression",
        'age_at_diagnosis',
        'time',
        'event'
    ]
    
    final_df = merged_df[cox_columns].copy()
    
    # 删除任何包含缺失值的行
    initial_rows = len(final_df)
    final_df.dropna(inplace=True)
    final_rows = len(final_df)
    print(f"  - 用于Cox分析的数据有 {initial_rows} 行。")
    print(f"  - 移除缺失值后，剩余 {final_rows} 行用于模型拟合。")

    if final_rows < 20: # 样本量太少，模型不可靠
        print("  - 错误: 清洗后剩余样本过少，无法进行可靠的Cox回归分析。")
        return

    # 4. 单变量Cox回归分析
    print("\n[4. 单变量Cox回归分析]")
    cph_univariate = CoxPHFitter()
    univariate_df = final_df[[f"{target_gene}_expression", 'time', 'event']]
    cph_univariate.fit(univariate_df, 'time', event_col='event')
    
    print("  - 单变量模型 (仅SMARCA2表达量) 结果:")
    cph_univariate.print_summary(model="univariate model", decimals=4)

    # 5. 多变量Cox回归分析
    print("\n[5. 多变量Cox回归分析]")
    cph_multivariate = CoxPHFitter()
    cph_multivariate.fit(final_df, 'time', event_col='event')
    
    print("  - 多变量模型 (校正年龄) 结果:")
    cph_multivariate.print_summary(model="multivariate model (age adjusted)", decimals=4)
    
    # 6. 可视化并保存结果
    print("\n[6. 生成并保存结果图表]")
    plt.figure(figsize=(10, 6))
    cph_multivariate.plot()
    plt.title(f'Cox Proportional Hazards Model - {target_gene}', fontsize=16, fontweight='bold')
    
    output_path = output_dir / f"{target_gene}_cox_regression_forest_plot"
    print(f"  - 正在保存森林图至 {output_path.name}.png")
    plt.savefig(f"{output_path}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_path}.pdf", bbox_inches='tight')
    plt.show()
    plt.close()
    
    print("\n=== Cox回归分析完成 ===")


def main():
    """
    Cox回归分析主流程。
    """
    parser = argparse.ArgumentParser(description="Perform Cox Proportional Hazards regression analysis.")
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

    # 创建输出目录
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    perform_cox_regression(
        expression_file=args.expression_file,
        clinical_file=args.clinical_file,
        target_gene=args.target_gene,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main() 