"""
直接绘制生存分析结果脚本
=============================

功能概述：
    从已有的批量生存分析结果文件中读取数据，
    直接绘制KM曲线和汇总图表，无需重新计算。

主要功能：
    1. 读取已有的分析结果文件
    2. 绘制前N个最显著基因的KM曲线
    3. 生成汇总统计图表
    4. 快速可视化展示

输入：
    - 生存分析结果CSV文件
    - 基因表达数据文件
    - 临床数据文件
    
输出：
    - KM曲线图（PNG/PDF格式）
    - 汇总统计图表

作者：BRM_Mitacs项目组
版本：1.0
"""

import argparse
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import logging
import sys
import re

# 忽略警告信息
warnings.filterwarnings('ignore')

# 定义项目根目录
PROJECT_ROOT = Path(__file__).resolve().parents[1]

def sanitize_filename(filename):
    """
    清理文件名中的非法字符
    
    参数：
        filename (str): 原始文件名
        
    返回：
        str: 清理后的安全文件名
    """
    # Windows文件名非法字符: < > : " | ? * / \
    illegal_chars = r'[<>:"|?*/\\]'
    # 将非法字符替换为下划线
    sanitized = re.sub(illegal_chars, '_', filename)
    # 移除连续的下划线
    sanitized = re.sub(r'_+', '_', sanitized)
    # 移除首尾的下划线和点号
    sanitized = sanitized.strip('_.')
    # 确保文件名不为空
    if not sanitized:
        sanitized = 'unnamed'
    return sanitized

def setup_logging(output_dir):
    """设置日志记录系统"""
    log_file = output_dir / "plot_survival_results_log.txt"
    
    logger = logging.getLogger('PlotSurvivalResults')
    logger.setLevel(logging.INFO)
    
    # 清除现有处理器
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 文件处理器
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.INFO)
    
    # 控制台处理器
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    
    # 格式化器
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

def load_data(expression_file, clinical_file, logger):
    """
    加载基因表达和临床数据
    
    参数：
        expression_file (Path): 基因表达数据文件路径
        clinical_file (Path): 临床数据文件路径
        logger: 日志记录器
        
    返回：
        tuple: (基因表达数据DataFrame, 生存数据DataFrame)
    """
    logger.info("正在加载基因表达和临床数据...")
    
    try:
        # 加载表达数据
        expr_df = pd.read_csv(expression_file, index_col=0)
        logger.info(f"✓ 表达数据维度: {expr_df.shape}")
        
        # 加载临床数据
        clin_df = pd.read_csv(clinical_file, sep='\t', index_col=0).T
        logger.info(f"✓ 临床数据维度: {clin_df.shape}")
        
        # 标准化样本ID
        expr_df.columns = expr_df.columns.str.upper().str.replace('.', '-').str.slice(0, 12)
        clin_df.index = clin_df.index.str.upper().str.slice(0, 12)
        
        # 去重
        expr_df = expr_df.loc[:, ~expr_df.columns.duplicated(keep='first')]
        clin_df = clin_df[~clin_df.index.duplicated(keep='first')]
        
        # 找到共同样本
        common_samples = list(set(expr_df.columns) & set(clin_df.index))
        logger.info(f"找到 {len(common_samples)} 个共同样本")
        
        # 准备生存数据
        survival_df = clin_df.loc[common_samples].copy()
        
        # 创建生存事件和时间列
        survival_df['event'] = pd.to_numeric(survival_df['vital_status'], errors='coerce')
        days_to_death = pd.to_numeric(survival_df['days_to_death'], errors='coerce')
        days_to_last_followup = pd.to_numeric(survival_df['days_to_last_followup'], errors='coerce')
        survival_df['time'] = np.where(survival_df['event'] == 1, days_to_death, days_to_last_followup)
        
        # 清理生存数据
        survival_df = survival_df.dropna(subset=['time', 'event'])
        survival_df = survival_df[survival_df['time'] > 0]
        
        # 筛选表达数据
        valid_samples = survival_df.index.tolist()
        expr_filtered = expr_df[valid_samples]
        
        logger.info(f"最终数据:")
        logger.info(f"  - 基因数量: {len(expr_filtered)}")
        logger.info(f"  - 样本数量: {len(expr_filtered.columns)}")
        logger.info(f"  - 事件数量 (死亡): {survival_df['event'].sum()}")
        logger.info(f"  - 删失数量 (存活): {(survival_df['event'] == 0).sum()}")
        
        return expr_filtered, survival_df
        
    except Exception as e:
        logger.error(f"✗ 数据加载失败: {e}")
        return None, None

def plot_km_curves_from_results(results_df, expr_df, survival_df, output_dir, top_n=10, logger=None):
    """
    根据分析结果绘制KM曲线
    
    参数：
        results_df (DataFrame): 生存分析结果
        expr_df (DataFrame): 基因表达数据
        survival_df (DataFrame): 生存数据
        output_dir (Path): 输出目录
        top_n (int): 绘制前N个最显著基因的KM曲线
        logger: 日志记录器
    """
    if logger:
        logger.info("="*80)
        logger.info(f"开始为前{top_n}个最显著基因绘制KM曲线")
        logger.info("="*80)
    
    # 选择前N个最显著的基因
    top_genes = results_df.head(top_n)
    
    # 设置图形样式
    plt.style.use('default')
    sns.set_palette("husl")
    
    for idx, (_, gene_info) in enumerate(top_genes.iterrows()):
        gene_name = gene_info['gene']
        
        # 确定使用哪个P值
        if 'cox_p_adjusted' in gene_info and pd.notna(gene_info['cox_p_adjusted']):
            p_value = gene_info['cox_p_adjusted']
            p_label = 'Cox P-value (Age Adjusted)'
        elif 'cox_p_unadjusted' in gene_info and pd.notna(gene_info['cox_p_unadjusted']):
            p_value = gene_info['cox_p_unadjusted']
            p_label = 'Cox P-value'
        elif 'logrank_p_value' in gene_info and pd.notna(gene_info['logrank_p_value']):
            p_value = gene_info['logrank_p_value']
            p_label = 'Log-rank P-value'
        else:
            p_value = gene_info.get('p_value', np.nan)
            p_label = 'P-value'
        
        if logger:
            logger.info(f"正在绘制基因 {gene_name} 的KM曲线 ({p_label}={p_value:.2e})")
        
        try:
            # 获取基因表达数据
            if gene_name not in expr_df.index:
                logger.warning(f"⚠️ 基因 {gene_name} 在表达数据中未找到，跳过")
                continue
                
            gene_data = expr_df.loc[gene_name]
            
            # 合并数据
            merged_df = pd.merge(
                pd.DataFrame({'expression': gene_data}),
                survival_df[['time', 'event']],
                left_index=True,
                right_index=True
            )
            
            # 分组
            low_quantile = merged_df['expression'].quantile(0.25)
            high_quantile = merged_df['expression'].quantile(0.75)
            
            high_group = merged_df[merged_df['expression'] >= high_quantile]
            low_group = merged_df[merged_df['expression'] <= low_quantile]
            
            if len(high_group) < 5 or len(low_group) < 5:
                logger.warning(f"⚠️ 基因 {gene_name} 分组样本数太少，跳过")
                continue
            
            # 绘制KM曲线
            plt.figure(figsize=(10, 8))
            
            kmf = KaplanMeierFitter()
            
            # 高表达组
            kmf.fit(
                durations=high_group['time'], 
                event_observed=high_group['event'], 
                label=f'High {gene_name} (n={len(high_group)})'
            )
            ax = kmf.plot(ci_show=False, linewidth=2.5)
            
            # 低表达组
            kmf.fit(
                durations=low_group['time'], 
                event_observed=low_group['event'], 
                label=f'Low {gene_name} (n={len(low_group)})'
            )
            kmf.plot(ax=ax, ci_show=False, linewidth=2.5)
            
            # 美化图表
            plt.title(f'Kaplan-Meier Survival Curve - {gene_name}', 
                     fontsize=16, fontweight='bold', pad=20)
            plt.xlabel('Time (days)', fontsize=14, fontweight='bold')
            plt.ylabel('Survival Probability', fontsize=14, fontweight='bold')
            plt.grid(True, linestyle='--', alpha=0.3)
            
            # 添加统计信息
            plt.text(0.02, 0.02, f'{p_label}: {p_value:.2e}', 
                    transform=ax.transAxes, fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
            
            # 添加排名信息
            plt.text(0.98, 0.98, f'Rank: #{idx+1}', 
                    transform=ax.transAxes, fontsize=12, fontweight='bold',
                    ha='right', va='top',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7))
            
            plt.legend(fontsize=12, frameon=True, fancybox=True, shadow=True)
            plt.tight_layout()
            
            # 保存图表
            safe_gene_name = sanitize_filename(gene_name)
            output_path = output_dir / f"KM_curve_{idx+1:02d}_{safe_gene_name}"
            plt.savefig(f"{output_path}.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{output_path}.pdf", bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            logger.error(f"✗ 绘制基因 {gene_name} 的KM曲线时出错: {e}")
            plt.close()
            continue
    
    if logger:
        logger.info(f"✓ 完成所有KM曲线绘制")

def create_summary_plots_from_results(results_df, output_dir, logger):
    """
    根据分析结果创建汇总图表
    
    参数：
        results_df (DataFrame): 生存分析结果
        output_dir (Path): 输出目录
        logger: 日志记录器
    """
    logger.info("正在创建汇总图表...")
    
    # 确定使用哪个P值列进行绘图
    if 'cox_p_adjusted' in results_df.columns and not results_df['cox_p_adjusted'].isna().all():
        p_col = 'cox_p_adjusted'
        p_title = 'Cox P-value (Age Adjusted)'
    elif 'cox_p_unadjusted' in results_df.columns and not results_df['cox_p_unadjusted'].isna().all():
        p_col = 'cox_p_unadjusted'
        p_title = 'Cox P-value (Unadjusted)'
    elif 'logrank_p_value' in results_df.columns:
        p_col = 'logrank_p_value'
        p_title = 'Log-rank P-value'
    else:
        p_col = 'p_value'  # 兼容旧格式
        p_title = 'P-value'
    
    # 1. P值分布直方图
    plt.figure(figsize=(16, 12))
    
    plt.subplot(2, 3, 1)
    p_values = results_df[p_col].dropna()
    plt.hist(p_values, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(x=0.05, color='red', linestyle='--', label='P=0.05')
    plt.axvline(x=0.01, color='orange', linestyle='--', label='P=0.01')
    plt.xlabel(p_title)
    plt.ylabel('基因数量')
    plt.title(f'{p_title}分布直方图')
    plt.legend()
    
    # 2. -log10(P值) vs 基因排名
    plt.subplot(2, 3, 2)
    ranks = range(1, len(p_values) + 1)
    plt.plot(ranks, -np.log10(p_values), 'o-', markersize=2, alpha=0.7)
    plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='P=0.05')
    plt.axhline(y=-np.log10(0.01), color='orange', linestyle='--', label='P=0.01')
    plt.xlabel('基因排名')
    plt.ylabel(f'-log10({p_title})')
    plt.title('显著性排名图')
    plt.legend()
    
    # 3. 样本数分布
    plt.subplot(2, 3, 3)
    total_samples = results_df['high_group_n'] + results_df['low_group_n']
    plt.hist(total_samples, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    plt.xlabel('分析样本数')
    plt.ylabel('基因数量')
    plt.title('样本数分布')
    
    # 4. 前20个最显著基因的条形图
    plt.subplot(2, 3, 4)
    top_20 = results_df.head(20)
    y_pos = range(len(top_20))
    plt.barh(y_pos, -np.log10(top_20[p_col]), color='coral', alpha=0.8)
    plt.yticks(y_pos, [gene[:15] + '...' if len(gene) > 15 else gene for gene in top_20['gene']])
    plt.xlabel(f'-log10({p_title})')
    plt.title('前20个最显著基因')
    plt.gca().invert_yaxis()
    
    # 5. 显著性分布饼图
    plt.subplot(2, 3, 5)
    sig_001 = len(results_df[results_df[p_col] < 0.01])
    sig_005 = len(results_df[(results_df[p_col] >= 0.01) & (results_df[p_col] < 0.05)])
    non_sig = len(results_df[results_df[p_col] >= 0.05])
    
    labels = ['P < 0.01', '0.01 ≤ P < 0.05', 'P ≥ 0.05']
    sizes = [sig_001, sig_005, non_sig]
    colors = ['red', 'orange', 'lightgray']
    
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    plt.title('显著性分布')
    
    # 6. 表达水平统计
    plt.subplot(2, 3, 6)
    if 'expression_mean' in results_df.columns:
        plt.scatter(results_df['expression_mean'], -np.log10(results_df[p_col]), alpha=0.5, s=10)
        plt.xlabel('平均表达水平')
        plt.ylabel(f'-log10({p_title})')
        plt.title('表达水平 vs 显著性')
    else:
        plt.text(0.5, 0.5, '表达水平数据不可用', ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('表达水平统计')
    
    plt.tight_layout()
    
    # 保存汇总图
    summary_path = output_dir / "survival_analysis_summary"
    plt.savefig(f"{summary_path}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{summary_path}.pdf", bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ 汇总图表已保存到 {summary_path.name}")

def main():
    """主函数"""
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="直接绘制生存分析结果 - 从已有结果文件生成图表")
    parser.add_argument(
        "--results_file",
        type=Path,
        default=PROJECT_ROOT / "results" / "batch_survival_analysis" / "all_genes_survival_analysis_results.csv",
        help="生存分析结果文件路径 (CSV格式)"
    )
    parser.add_argument(
        "--expression_file",
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "working.csv",
        help="基因表达数据文件路径 (CSV格式)"
    )
    parser.add_argument(
        "--clinical_file", 
        type=Path,
        default=PROJECT_ROOT / "data" / "processed" / "OV.clin.merged.picked.txt",
        help="临床数据文件路径"
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "batch_survival_analysis",
        help="输出目录"
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=10,
        help="绘制前N个最显著基因的KM曲线 (默认: 10)"
    )
    
    args = parser.parse_args()
    
    # 创建输出目录
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # 设置日志
    logger = setup_logging(args.output_dir)
    
    # 记录开始时间
    start_time = datetime.now()
    logger.info("="*80)
    logger.info("直接绘制生存分析结果开始")
    logger.info("="*80)
    logger.info(f"绘制开始时间: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"配置参数:")
    logger.info(f"  - 结果文件: {args.results_file.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 表达数据文件: {args.expression_file.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 临床数据文件: {args.clinical_file.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 输出目录: {args.output_dir.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 绘制KM曲线数: {args.top_n}")
    
    try:
        # 1. 加载分析结果
        logger.info("="*80)
        logger.info("加载分析结果")
        logger.info("="*80)
        
        if not args.results_file.exists():
            logger.error(f"✗ 结果文件不存在: {args.results_file}")
            return
        
        results_df = pd.read_csv(args.results_file)
        logger.info(f"✓ 成功加载分析结果: {len(results_df)} 个基因")
        
        # 2. 加载原始数据（用于绘制KM曲线）
        logger.info("="*80)
        logger.info("加载原始数据")
        logger.info("="*80)
        
        expr_df, survival_df = load_data(args.expression_file, args.clinical_file, logger)
        if expr_df is None or survival_df is None:
            logger.error("✗ 原始数据加载失败")
            return
        
        # 3. 绘制KM曲线
        logger.info("="*80)
        logger.info("绘制KM曲线")
        logger.info("="*80)
        
        plot_km_curves_from_results(
            results_df, expr_df, survival_df, 
            args.output_dir, args.top_n, logger
        )
        
        # 4. 创建汇总图表
        logger.info("="*80)
        logger.info("创建汇总图表")
        logger.info("="*80)
        
        create_summary_plots_from_results(results_df, args.output_dir, logger)
        
        # 记录结束时间
        end_time = datetime.now()
        total_duration = (end_time - start_time).total_seconds()
        
        logger.info("="*80)
        logger.info("直接绘制生存分析结果完成")
        logger.info("="*80)
        logger.info(f"绘制结束时间: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"总耗时: {total_duration:.1f} 秒")
        
        # 输出结果总结
        logger.info(f"\n结果总结:")
        logger.info(f"  - 分析基因数: {len(results_df)}")
        
        # 确定使用哪个P值列
        if 'cox_p_adjusted' in results_df.columns and not results_df['cox_p_adjusted'].isna().all():
            p_col = 'cox_p_adjusted'
        elif 'cox_p_unadjusted' in results_df.columns and not results_df['cox_p_unadjusted'].isna().all():
            p_col = 'cox_p_unadjusted'
        elif 'logrank_p_value' in results_df.columns:
            p_col = 'logrank_p_value'
        else:
            p_col = 'p_value'
        
        significant_genes = results_df[results_df[p_col] < 0.05]
        logger.info(f"  - 显著基因数(P<0.05): {len(significant_genes)}")
        logger.info(f"  - 生成KM曲线数: {min(args.top_n, len(results_df))}")
        logger.info(f"  - 最小P值: {results_df[p_col].min():.2e}")
        
        logger.info(f"\n所有图表已保存到: {args.output_dir.relative_to(PROJECT_ROOT)}")
        
    except Exception as e:
        logger.error(f"✗ 绘制过程中发生错误: {e}")
        import traceback
        logger.error(f"错误详情: {traceback.format_exc()}")
        return

if __name__ == "__main__":
    main() 