"""
批量生存分析脚本
================

功能概述：
    对所有基因进行批量生存分析，识别与患者生存预后最相关的基因。
    使用Log-rank检验评估每个基因的生存分析显著性，
    并为最显著的基因生成Kaplan-Meier生存曲线图。

主要功能：
    1. 批量处理所有基因的生存分析
    2. 使用Log-rank检验评估统计显著性
    3. 筛选最显著的基因（Top N）
    4. 生成详细的生存分析报告
    5. 为重要基因创建高质量的KM曲线图

输入：
    - 基因表达数据文件（CSV格式）
    - 临床数据文件（制表符分隔）
    
输出：
    - 所有基因的生存分析结果表格
    - 最显著基因的KM曲线图
    - 详细的分析日志和统计报告

作者：BRM_Mitacs项目组
版本：1.0
"""

import argparse
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import logging
import sys
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import re

# 忽略警告信息
warnings.filterwarnings('ignore')

# 导入Cox回归相关库
from lifelines import CoxPHFitter

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
    """
    设置日志记录系统
    
    参数：
        output_dir (Path): 日志文件保存目录
    """
    log_file = output_dir / "batch_survival_analysis_log.txt"
    
    logger = logging.getLogger('BatchSurvivalAnalysis')
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

def prepare_survival_data(expression_file, clinical_file, logger):
    """
    准备生存分析数据
    
    参数：
        expression_file (Path): 基因表达数据文件路径
        clinical_file (Path): 临床数据文件路径
        logger: 日志记录器
        
    返回：
        tuple: (基因表达数据DataFrame, 合并的生存数据DataFrame)
    """
    logger.info("="*80)
    logger.info("开始数据加载和预处理")
    logger.info("="*80)
    
    # 1. 加载数据
    logger.info("正在加载基因表达和临床数据...")
    try:
        expr_df = pd.read_csv(expression_file, index_col=0)
        clin_df = pd.read_csv(clinical_file, sep='\t', index_col=0).T
        logger.info(f"✓ 原始表达数据维度: {expr_df.shape}")
        logger.info(f"✓ 原始临床数据维度: {clin_df.shape}")
    except FileNotFoundError as e:
        logger.error(f"✗ 数据加载错误: {e}")
        return None, None
    
    # 2. 标准化样本ID
    logger.info("正在标准化样本ID...")
    expr_df.columns = expr_df.columns.str.upper().str.replace('.', '-').str.slice(0, 12)
    clin_df.index = clin_df.index.str.upper().str.slice(0, 12)
    
    # 检查并显示重复样本
    original_expr_samples = len(expr_df.columns)
    original_clin_samples = len(clin_df)
    
    # 找出表达数据中的重复样本
    expr_duplicated = expr_df.columns.duplicated(keep=False)  # 标记所有重复项
    if expr_duplicated.any():
        duplicated_expr_samples = expr_df.columns[expr_duplicated].unique()
        logger.info(f"发现表达数据中的重复样本 ({len(duplicated_expr_samples)} 个):")
        for sample in duplicated_expr_samples:
            count = (expr_df.columns == sample).sum()
            logger.info(f"  - {sample}: 出现 {count} 次")
    else:
        logger.info("表达数据中未发现重复样本")
    
    # 找出临床数据中的重复样本
    clin_duplicated = clin_df.index.duplicated(keep=False)  # 标记所有重复项
    if clin_duplicated.any():
        duplicated_clin_samples = clin_df.index[clin_duplicated].unique()
        logger.info(f"发现临床数据中的重复样本 ({len(duplicated_clin_samples)} 个):")
        for sample in duplicated_clin_samples:
            count = (clin_df.index == sample).sum()
            logger.info(f"  - {sample}: 出现 {count} 次")
    else:
        logger.info("临床数据中未发现重复样本")
    
    # 去重（保留第一个出现的样本）
    expr_df = expr_df.loc[:, ~expr_df.columns.duplicated(keep='first')]
    clin_df = clin_df[~clin_df.index.duplicated(keep='first')]
    
    logger.info(f"表达数据去重后: {len(expr_df.columns)} / {original_expr_samples} 个唯一样本")
    logger.info(f"临床数据去重后: {len(clin_df)} / {original_clin_samples} 个唯一样本")
    
    # 3. 找到共同样本
    common_samples = list(set(expr_df.columns) & set(clin_df.index))
    logger.info(f"找到 {len(common_samples)} 个共同样本")
    
    if len(common_samples) < 50:
        logger.error(f"✗ 共同样本数量太少 ({len(common_samples)})，无法进行可靠的生存分析")
        return None, None
    
    # 4. 准备生存数据
    logger.info("正在处理生存数据...")
    survival_df = clin_df.loc[common_samples].copy()
    
    # 创建生存事件和时间列
    survival_df['event'] = pd.to_numeric(survival_df['vital_status'], errors='coerce')
    days_to_death = pd.to_numeric(survival_df['days_to_death'], errors='coerce')
    days_to_last_followup = pd.to_numeric(survival_df['days_to_last_followup'], errors='coerce')
    survival_df['time'] = np.where(survival_df['event'] == 1, days_to_death, days_to_last_followup)
    
    # 清理生存数据
    original_survival_samples = len(survival_df)
    survival_df = survival_df.dropna(subset=['time', 'event'])
    survival_df = survival_df[survival_df['time'] > 0]  # 移除时间为0或负数的样本
    
    logger.info(f"生存数据清理: {len(survival_df)} / {original_survival_samples} 个样本有完整生存数据")
    
    # 5. 筛选表达数据
    valid_samples = survival_df.index.tolist()
    expr_filtered = expr_df[valid_samples]
    
    # 6. 处理年龄信息
    logger.info("正在处理年龄信息...")
    
    # 检查年龄数据
    if 'years_to_birth' in survival_df.columns:
        survival_df['age'] = pd.to_numeric(survival_df['years_to_birth'], errors='coerce')
        
        # 检查年龄数据完整性
        age_missing = survival_df['age'].isna().sum()
        age_valid = len(survival_df) - age_missing
        
        logger.info(f"年龄数据统计:")
        logger.info(f"  - 有效年龄数据: {age_valid} / {len(survival_df)} 个样本")
        logger.info(f"  - 缺失年龄数据: {age_missing} 个样本")
        logger.info(f"  - 年龄范围: [{survival_df['age'].min():.0f}, {survival_df['age'].max():.0f}] 岁")
        logger.info(f"  - 平均年龄: {survival_df['age'].mean():.1f} ± {survival_df['age'].std():.1f} 岁")
        logger.info(f"  - 中位年龄: {survival_df['age'].median():.0f} 岁")
        
        # 年龄分组分析
        young_threshold = survival_df['age'].quantile(0.33)  # 年轻组（前1/3）
        old_threshold = survival_df['age'].quantile(0.67)    # 年老组（后1/3）
        
        survival_df['age_group'] = 'middle'
        survival_df.loc[survival_df['age'] <= young_threshold, 'age_group'] = 'young'
        survival_df.loc[survival_df['age'] >= old_threshold, 'age_group'] = 'old'
        
        age_group_counts = survival_df['age_group'].value_counts()
        logger.info(f"  - 年龄分组: 年轻组≤{young_threshold:.0f}岁({age_group_counts.get('young', 0)}人), "
                   f"中年组({age_group_counts.get('middle', 0)}人), "
                   f"年老组≥{old_threshold:.0f}岁({age_group_counts.get('old', 0)}人)")
        
        # 检查年龄与生存的关系
        if age_valid >= 50:  # 只有足够样本时才进行分析
            from lifelines.statistics import logrank_test
            young_group = survival_df[survival_df['age_group'] == 'young']
            old_group = survival_df[survival_df['age_group'] == 'old']
            
            if len(young_group) >= 10 and len(old_group) >= 10:
                age_survival_test = logrank_test(
                    young_group['time'], old_group['time'],
                    young_group['event'], old_group['event']
                )
                logger.info(f"  - 年龄与生存关联性检验: P = {age_survival_test.p_value:.3f}")
                if age_survival_test.p_value < 0.05:
                    logger.warning(f"  ⚠️ 年龄显著影响生存预后 (P < 0.05)，建议在分析中调整年龄因素")
                else:
                    logger.info(f"  ✓ 年龄对生存预后影响不显著 (P ≥ 0.05)")
        
        # 移除年龄缺失的样本
        if age_missing > 0:
            logger.info(f"移除 {age_missing} 个年龄数据缺失的样本")
            survival_df = survival_df.dropna(subset=['age'])
            valid_samples = survival_df.index.tolist()
            expr_filtered = expr_df[valid_samples]
    else:
        logger.warning("临床数据中未找到年龄信息 (years_to_birth)，将无法调整年龄因素")
        survival_df['age'] = None
    
    logger.info(f"最终用于分析的数据:")
    logger.info(f"  - 基因数量: {len(expr_filtered)}")
    logger.info(f"  - 样本数量: {len(expr_filtered.columns)}")
    logger.info(f"  - 事件数量 (死亡): {survival_df['event'].sum()}")
    logger.info(f"  - 删失数量 (存活): {(survival_df['event'] == 0).sum()}")
    
    # 7. 数据质量检查和预览保存
    logger.info("="*80)
    logger.info("数据质量检查")
    logger.info("="*80)
    
    # 检查表达数据分布
    logger.info("表达数据统计:")
    logger.info(f"  - 表达值范围: [{expr_filtered.values.min():.4f}, {expr_filtered.values.max():.4f}]")
    logger.info(f"  - 表达值均值: {expr_filtered.values.mean():.4f}")
    logger.info(f"  - 表达值标准差: {expr_filtered.values.std():.4f}")
    logger.info(f"  - 零值比例: {(expr_filtered.values == 0).sum() / expr_filtered.size * 100:.2f}%")
    
    # 检查生存时间分布
    logger.info("生存时间统计:")
    logger.info(f"  - 时间范围: [{survival_df['time'].min():.0f}, {survival_df['time'].max():.0f}] 天")
    logger.info(f"  - 中位随访时间: {survival_df['time'].median():.0f} 天")
    logger.info(f"  - 平均随访时间: {survival_df['time'].mean():.0f} 天")
    
    # 保存数据预览
    output_dir = PROJECT_ROOT / "results" / "batch_survival_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("="*80)
    logger.info("保存数据预览文件")
    logger.info("="*80)
    
    # 保存生存数据预览
    survival_preview = survival_df[['time', 'event', 'vital_status', 'days_to_death', 'days_to_last_followup']].copy()
    survival_preview_file = output_dir / "survival_data_preview.csv"
    survival_preview.to_csv(survival_preview_file)
    logger.info(f"✓ 生存数据预览已保存: {survival_preview_file.name}")
    
    # 保存前10个基因的表达数据预览
    expr_preview = expr_filtered.head(10).T  # 转置，样本为行，基因为列
    expr_preview_file = output_dir / "expression_data_preview_top10_genes.csv"
    expr_preview.to_csv(expr_preview_file)
    logger.info(f"✓ 表达数据预览(前10个基因)已保存: {expr_preview_file.name}")
    
    # 保存数据维度信息
    data_info = {
        'metric': [
            'total_genes', 'total_samples', 'death_events', 'censored_events',
            'expr_min', 'expr_max', 'expr_mean', 'expr_std', 'zero_proportion_pct',
            'time_min_days', 'time_max_days', 'time_median_days', 'time_mean_days'
        ],
        'value': [
            len(expr_filtered), len(expr_filtered.columns), survival_df['event'].sum(), 
            (survival_df['event'] == 0).sum(), expr_filtered.values.min(), expr_filtered.values.max(),
            expr_filtered.values.mean(), expr_filtered.values.std(), 
            (expr_filtered.values == 0).sum() / expr_filtered.size * 100,
            survival_df['time'].min(), survival_df['time'].max(),
            survival_df['time'].median(), survival_df['time'].mean()
        ]
    }
    data_info_df = pd.DataFrame(data_info)
    data_info_file = output_dir / "data_summary_statistics.csv"
    data_info_df.to_csv(data_info_file, index=False)
    logger.info(f"✓ 数据统计信息已保存: {data_info_file.name}")
    
    # 检查是否有表达数据为常数的基因
    constant_genes = []
    for gene in expr_filtered.index:
        if expr_filtered.loc[gene].std() == 0:
            constant_genes.append(gene)
    
    if constant_genes:
        logger.warning(f"发现 {len(constant_genes)} 个表达值为常数的基因，将被排除")
        expr_filtered = expr_filtered.drop(constant_genes)
        logger.info(f"排除常数基因后，剩余 {len(expr_filtered)} 个基因")
    
    logger.info("="*80)
    logger.info("数据预处理完成，可以开始分析")
    logger.info("="*80)
    
    return expr_filtered, survival_df

def analyze_single_gene_with_age_adjustment(gene_data, survival_df, gene_name):
    """
    对单个基因进行调整年龄的生存分析
    
    参数：
        gene_data (Series): 单个基因的表达数据
        survival_df (DataFrame): 生存数据（包含年龄信息）
        gene_name (str): 基因名称
        
    返回：
        dict: 分析结果字典
    """
    try:
        # 合并基因表达和生存数据
        merged_df = pd.merge(
            pd.DataFrame({'expression': gene_data}),
            survival_df[['time', 'event', 'age']],
            left_index=True,
            right_index=True
        )
        
        if len(merged_df) < 20:  # 样本数太少
            return None
        
        # 检查是否有年龄数据
        has_age = not merged_df['age'].isna().all()
        
        # 根据分位数分组（低表达组：最低40%，高表达组：最高40%）
        low_quantile = merged_df['expression'].quantile(0.40)
        high_quantile = merged_df['expression'].quantile(0.60)
        
        high_group = merged_df[merged_df['expression'] >= high_quantile]
        low_group = merged_df[merged_df['expression'] <= low_quantile]
        
        if len(high_group) < 5 or len(low_group) < 5:
            return None
        
        # 1. Log-rank检验（未调整年龄）
        try:
            logrank_results = logrank_test(
                durations_A=high_group['time'], 
                event_observed_A=high_group['event'],
                durations_B=low_group['time'], 
                event_observed_B=low_group['event']
            )
            logrank_p = logrank_results.p_value
        except Exception:
            logrank_p = np.nan
        
        # 2. Cox回归分析（调整年龄）
        cox_p_unadjusted = np.nan
        cox_p_adjusted = np.nan
        cox_hr_unadjusted = np.nan
        cox_hr_adjusted = np.nan
        cox_hr_ci_lower_unadjusted = np.nan
        cox_hr_ci_upper_unadjusted = np.nan
        cox_hr_ci_lower_adjusted = np.nan
        cox_hr_ci_upper_adjusted = np.nan
        
        try:
            # 准备Cox回归数据（连续变量）
            cox_df = merged_df.copy()
            cox_df = cox_df.dropna()
            
            if len(cox_df) >= 20:
                # 单变量Cox回归（仅基因表达）
                try:
                    cph_univariate = CoxPHFitter()
                    cph_univariate.fit(cox_df[['time', 'event', 'expression']], 
                                     duration_col='time', event_col='event')
                    
                    cox_p_unadjusted = cph_univariate.summary.loc['expression', 'p']
                    cox_hr_unadjusted = cph_univariate.summary.loc['expression', 'exp(coef)']
                    cox_hr_ci_lower_unadjusted = cph_univariate.summary.loc['expression', 'exp(coef) lower 95%']
                    cox_hr_ci_upper_unadjusted = cph_univariate.summary.loc['expression', 'exp(coef) upper 95%']
                except Exception:
                    pass
                
                # 多变量Cox回归（基因表达 + 年龄）
                if has_age and not cox_df['age'].isna().all():
                    try:
                        cph_multivariate = CoxPHFitter()
                        cph_multivariate.fit(cox_df[['time', 'event', 'expression', 'age']], 
                                           duration_col='time', event_col='event')
                        
                        cox_p_adjusted = cph_multivariate.summary.loc['expression', 'p']
                        cox_hr_adjusted = cph_multivariate.summary.loc['expression', 'exp(coef)']
                        cox_hr_ci_lower_adjusted = cph_multivariate.summary.loc['expression', 'exp(coef) lower 95%']
                        cox_hr_ci_upper_adjusted = cph_multivariate.summary.loc['expression', 'exp(coef) upper 95%']
                    except Exception:
                        pass
        except Exception:
            pass
        
        # 计算中位生存时间
        kmf = KaplanMeierFitter()
        
        try:
            kmf.fit(high_group['time'], high_group['event'])
            median_survival_high = kmf.median_survival_time_
        except Exception:
            median_survival_high = np.nan
        
        try:
            kmf.fit(low_group['time'], low_group['event'])
            median_survival_low = kmf.median_survival_time_
        except Exception:
            median_survival_low = np.nan
        
        # 年龄统计
        if has_age:
            age_mean_high = high_group['age'].mean()
            age_mean_low = low_group['age'].mean()
            age_p_value = np.nan
            
            # 检验高低表达组年龄差异
            try:
                from scipy.stats import ttest_ind
                _, age_p_value = ttest_ind(high_group['age'].dropna(), 
                                        low_group['age'].dropna())
            except Exception:
                pass
        else:
            age_mean_high = np.nan
            age_mean_low = np.nan
            age_p_value = np.nan
        
        return {
            'gene': gene_name,
            # Log-rank检验结果
            'logrank_p_value': logrank_p,
            'logrank_test_statistic': logrank_results.test_statistic if 'logrank_results' in locals() else np.nan,
            # Cox回归结果（未调整）
            'cox_p_unadjusted': cox_p_unadjusted,
            'cox_hr_unadjusted': cox_hr_unadjusted,
            'cox_hr_ci_lower_unadjusted': cox_hr_ci_lower_unadjusted,
            'cox_hr_ci_upper_unadjusted': cox_hr_ci_upper_unadjusted,
            # Cox回归结果（年龄调整）
            'cox_p_adjusted': cox_p_adjusted,
            'cox_hr_adjusted': cox_hr_adjusted,
            'cox_hr_ci_lower_adjusted': cox_hr_ci_lower_adjusted,
            'cox_hr_ci_upper_adjusted': cox_hr_ci_upper_adjusted,
            # 基本统计信息
            'high_group_n': len(high_group),
            'low_group_n': len(low_group),
            'high_group_events': high_group['event'].sum(),
            'low_group_events': low_group['event'].sum(),
            'median_survival_high': median_survival_high,
            'median_survival_low': median_survival_low,
            'expression_mean': merged_df['expression'].mean(),
            'expression_std': merged_df['expression'].std(),
            # 年龄相关统计
            'age_mean_high': age_mean_high,
            'age_mean_low': age_mean_low,
            'age_difference_p_value': age_p_value,
            'has_age_data': has_age
        }
        
    except Exception as e:
        return None

def analyze_single_gene(gene_data, survival_df, gene_name):
    """
    对单个基因进行生存分析（保持向后兼容性）
    """
    # 如果有年龄数据，使用调整年龄的分析
    if 'age' in survival_df.columns and not survival_df['age'].isna().all():
        return analyze_single_gene_with_age_adjustment(gene_data, survival_df, gene_name)
    
    # 否则使用原始的分析方法
    try:
        # 合并基因表达和生存数据
        merged_df = pd.merge(
            pd.DataFrame({'expression': gene_data}),
            survival_df[['time', 'event']],
            left_index=True,
            right_index=True
        )
        
        if len(merged_df) < 20:  # 样本数太少
            return None
        
        # 根据分位数分组（低表达组：最低40%，高表达组：最高40%）
        low_quantile = merged_df['expression'].quantile(0.40)
        high_quantile = merged_df['expression'].quantile(0.60)
        
        high_group = merged_df[merged_df['expression'] >= high_quantile]
        low_group = merged_df[merged_df['expression'] <= low_quantile]
        
        if len(high_group) < 5 or len(low_group) < 5:
            return None
        
        # 进行Log-rank检验
        try:
            results = logrank_test(
                durations_A=high_group['time'], 
                event_observed_A=high_group['event'],
                durations_B=low_group['time'], 
                event_observed_B=low_group['event']
            )
            
            # 计算中位生存时间
            kmf = KaplanMeierFitter()
            
            kmf.fit(high_group['time'], high_group['event'])
            median_survival_high = kmf.median_survival_time_
            
            kmf.fit(low_group['time'], low_group['event'])
            median_survival_low = kmf.median_survival_time_
            
            return {
                'gene': gene_name,
                'p_value': results.p_value,
                'test_statistic': results.test_statistic,
                'high_group_n': len(high_group),
                'low_group_n': len(low_group),
                'high_group_events': high_group['event'].sum(),
                'low_group_events': low_group['event'].sum(),
                'median_survival_high': median_survival_high,
                'median_survival_low': median_survival_low,
                'expression_mean': merged_df['expression'].mean(),
                'expression_std': merged_df['expression'].std()
            }
            
        except Exception as e:
            return None
            
    except Exception as e:
        return None

def batch_survival_analysis(expr_df, survival_df, logger, n_cores=None):
    """
    批量进行生存分析
    
    参数：
        expr_df (DataFrame): 基因表达数据
        survival_df (DataFrame): 生存数据
        logger: 日志记录器
        n_cores (int): 使用的CPU核心数
        
    返回：
        DataFrame: 所有基因的生存分析结果
    """
    logger.info("="*80)
    logger.info("开始批量生存分析")
    logger.info("="*80)
    
    # 设置并行处理 - 允许使用CPU总数-1个核心
    if n_cores is None:
        n_cores = mp.cpu_count() - 1  # 使用CPU总数-1个核心
    else:
        n_cores = min(n_cores, mp.cpu_count() - 1)  # 限制不超过CPU总数-1
    
    logger.info(f"系统可用CPU核心数: {mp.cpu_count()}")
    logger.info(f"实际使用CPU核心数: {n_cores}")
    logger.info(f"总共需要分析基因数: {len(expr_df)}")
    
    # 估算分析时间
    estimated_time_per_gene = 0.1  # 每个基因大约0.1秒
    total_estimated_time = len(expr_df) * estimated_time_per_gene / n_cores
    logger.info(f"预估总分析时间: {total_estimated_time/60:.1f} 分钟")
    
    # 准备分析函数
    analyze_func = partial(analyze_single_gene, survival_df=survival_df)
    
    # 准备基因数据
    logger.info("正在准备基因数据...")
    gene_list = []
    failed_genes = []
    
    for i, gene_name in enumerate(expr_df.index):
        try:
            gene_data = expr_df.loc[gene_name]
            gene_data.name = gene_name
            
            # 基本数据质量检查
            if gene_data.isna().all():
                failed_genes.append((gene_name, "所有值为NaN"))
                continue
            if gene_data.std() == 0:
                failed_genes.append((gene_name, "表达值为常数"))
                continue
                
            gene_list.append((gene_data, survival_df, gene_name))
            
        except Exception as e:
            failed_genes.append((gene_name, f"数据准备错误: {str(e)}"))
    
    logger.info(f"数据准备完成:")
    logger.info(f"  - 可分析基因数: {len(gene_list)}")
    logger.info(f"  - 跳过基因数: {len(failed_genes)}")
    
    if failed_genes:
        logger.info(f"跳过的基因示例 (前5个): {[f'{gene}({reason})' for gene, reason in failed_genes[:5]]}")
    
    if len(gene_list) == 0:
        logger.error("✗ 没有基因可以进行分析")
        return None
    
    # 使用多进程进行分析
    logger.info("="*50)
    logger.info("开始并行分析进程")
    logger.info("="*50)
    
    results = []
    successful_analyses = 0
    failed_analyses = 0
    
    # 分批处理以提供更好的进度反馈
    batch_size = max(100, len(gene_list) // 20)  # 分成大约20个批次
    total_batches = (len(gene_list) + batch_size - 1) // batch_size
    
    logger.info(f"分析配置:")
    logger.info(f"  - 批次大小: {batch_size} 个基因/批次")
    logger.info(f"  - 总批次数: {total_batches}")
    logger.info(f"  - 并行进程数: {n_cores}")
    
    start_time = datetime.now()
    
    with mp.Pool(processes=n_cores) as pool:
        for batch_idx in range(total_batches):
            batch_start = batch_idx * batch_size
            batch_end = min((batch_idx + 1) * batch_size, len(gene_list))
            batch_genes = gene_list[batch_start:batch_end]
            
            logger.info(f"\n批次 {batch_idx + 1}/{total_batches}: 分析基因 {batch_start + 1}-{batch_end}")
            
            # 分析当前批次
            batch_results = list(tqdm(
                pool.starmap(analyze_single_gene, batch_genes),
                total=len(batch_genes),
                desc=f"批次 {batch_idx + 1}/{total_batches}",
                ncols=100,
                leave=False
            ))
            
            # 统计当前批次结果
            batch_successful = sum(1 for r in batch_results if r is not None)
            batch_failed = len(batch_results) - batch_successful
            successful_analyses += batch_successful
            failed_analyses += batch_failed
            
            results.extend(batch_results)
            
            # 批次完成统计
            elapsed_time = (datetime.now() - start_time).total_seconds()
            progress_pct = (batch_end / len(gene_list)) * 100
            estimated_remaining = elapsed_time / progress_pct * (100 - progress_pct) if progress_pct > 0 else 0
            
            logger.info(f"  ✓ 批次完成: {batch_successful}/{len(batch_genes)} 个基因成功分析")
            logger.info(f"  总进度: {progress_pct:.1f}% ({successful_analyses}/{len(gene_list)})")
            logger.info(f"  已用时间: {elapsed_time/60:.1f} 分钟")
            logger.info(f"  预计剩余: {estimated_remaining/60:.1f} 分钟")
    
    total_time = (datetime.now() - start_time).total_seconds()
    
    logger.info("="*50)
    logger.info("并行分析完成")
    logger.info("="*50)
    logger.info(f"分析统计:")
    logger.info(f"  - 成功分析: {successful_analyses} 个基因")
    logger.info(f"  - 分析失败: {failed_analyses} 个基因")
    logger.info(f"  - 成功率: {successful_analyses/(successful_analyses+failed_analyses)*100:.1f}%")
    logger.info(f"  - 总耗时: {total_time/60:.1f} 分钟")
    logger.info(f"  - 平均速度: {successful_analyses/total_time:.1f} 基因/秒")
    
    # 过滤掉None结果
    valid_results = [r for r in results if r is not None]
    
    if len(valid_results) == 0:
        logger.error("✗ 没有基因通过生存分析")
        return None
    
    # 转换为DataFrame
    logger.info("正在处理分析结果...")
    results_df = pd.DataFrame(valid_results)
    
    # 添加多重检验校正
    logger.info("正在进行多重检验校正...")
    from statsmodels.stats.multitest import fdrcorrection
    
    # 检查结果格式，决定使用哪个P值进行校正
    if 'cox_p_adjusted' in results_df.columns and not results_df['cox_p_adjusted'].isna().all():
        # 使用年龄调整后的Cox回归P值
        valid_p_values = results_df['cox_p_adjusted'].dropna()
        if len(valid_p_values) > 0:
            _, results_df['cox_p_adjusted_fdr'] = fdrcorrection(results_df['cox_p_adjusted'].fillna(1), alpha=0.05)
            primary_p_col = 'cox_p_adjusted'
            logger.info("使用年龄调整后的Cox回归P值进行多重检验校正")
        else:
            primary_p_col = 'logrank_p_value'
    elif 'cox_p_unadjusted' in results_df.columns and not results_df['cox_p_unadjusted'].isna().all():
        # 使用未调整的Cox回归P值
        _, results_df['cox_p_unadjusted_fdr'] = fdrcorrection(results_df['cox_p_unadjusted'].fillna(1), alpha=0.05)
        primary_p_col = 'cox_p_unadjusted'
        logger.info("使用未调整的Cox回归P值进行多重检验校正")
    elif 'logrank_p_value' in results_df.columns:
        # 使用Log-rank检验P值
        _, results_df['logrank_p_fdr'] = fdrcorrection(results_df['logrank_p_value'].fillna(1), alpha=0.05)
        primary_p_col = 'logrank_p_value'
        logger.info("使用Log-rank检验P值进行多重检验校正")
    else:
        # 兼容旧格式
        _, results_df['p_value_adj'] = fdrcorrection(results_df['p_value'].fillna(1), alpha=0.05)
        primary_p_col = 'p_value'
        logger.info("使用传统P值进行多重检验校正")
    
    # 排序
    results_df = results_df.sort_values(primary_p_col).reset_index(drop=True)
    
    # 输出统计信息
    if primary_p_col == 'cox_p_adjusted':
        significant_genes = results_df[results_df['cox_p_adjusted'] < 0.05]
        highly_significant_genes = results_df[results_df['cox_p_adjusted'] < 0.01]
        adj_significant_genes = results_df[results_df['cox_p_adjusted_fdr'] < 0.05]
        min_p_value = results_df['cox_p_adjusted'].min()
        median_p_value = results_df['cox_p_adjusted'].median()
    elif primary_p_col == 'cox_p_unadjusted':
        significant_genes = results_df[results_df['cox_p_unadjusted'] < 0.05]
        highly_significant_genes = results_df[results_df['cox_p_unadjusted'] < 0.01]
        adj_significant_genes = results_df[results_df['cox_p_unadjusted_fdr'] < 0.05]
        min_p_value = results_df['cox_p_unadjusted'].min()
        median_p_value = results_df['cox_p_unadjusted'].median()
    elif primary_p_col == 'logrank_p_value':
        significant_genes = results_df[results_df['logrank_p_value'] < 0.05]
        highly_significant_genes = results_df[results_df['logrank_p_value'] < 0.01]
        adj_significant_genes = results_df[results_df['logrank_p_fdr'] < 0.05]
        min_p_value = results_df['logrank_p_value'].min()
        median_p_value = results_df['logrank_p_value'].median()
    else:
        # 兼容旧格式
        significant_genes = results_df[results_df['p_value'] < 0.05]
        highly_significant_genes = results_df[results_df['p_value'] < 0.01]
        adj_significant_genes = results_df[results_df['p_value_adj'] < 0.05]
        min_p_value = results_df['p_value'].min()
        median_p_value = results_df['p_value'].median()
    
    logger.info("="*80)
    logger.info("生存分析统计结果")
    logger.info("="*80)
    logger.info(f"结果概览:")
    logger.info(f"  - 总分析基因数: {len(valid_results)}")
    logger.info(f"  - P < 0.05 的基因数: {len(significant_genes)} ({len(significant_genes)/len(valid_results)*100:.1f}%)")
    logger.info(f"  - P < 0.01 的基因数: {len(highly_significant_genes)} ({len(highly_significant_genes)/len(valid_results)*100:.1f}%)")
    logger.info(f"  - 校正后P < 0.05 的基因数: {len(adj_significant_genes)} ({len(adj_significant_genes)/len(valid_results)*100:.1f}%)")
    logger.info(f"  - 最小P值: {min_p_value:.2e}")
    logger.info(f"  - P值中位数: {median_p_value:.3f}")
    
    # 年龄调整效果分析
    if 'cox_p_adjusted' in results_df.columns and 'cox_p_unadjusted' in results_df.columns:
        # 比较调整前后的显著基因数
        unadj_sig = len(results_df[results_df['cox_p_unadjusted'] < 0.05])
        adj_sig = len(results_df[results_df['cox_p_adjusted'] < 0.05])
        
        logger.info(f"\n年龄调整效果:")
        logger.info(f"  - 未调整年龄时显著基因数: {unadj_sig}")
        logger.info(f"  - 调整年龄后显著基因数: {adj_sig}")
        logger.info(f"  - 变化: {adj_sig - unadj_sig:+d} 个基因")
        
        if adj_sig < unadj_sig:
            logger.info(f"  💡 年龄调整减少了 {unadj_sig - adj_sig} 个假阳性结果")
        elif adj_sig > unadj_sig:
            logger.info(f"  💡 年龄调整发现了 {adj_sig - unadj_sig} 个新的显著基因")
        else:
            logger.info(f"  💡 年龄调整对结果影响较小")
    
    # 显示前10个最显著的基因
    logger.info(f"\n前10个最显著的基因:")
    for i, (_, gene_info) in enumerate(results_df.head(10).iterrows()):
        if primary_p_col == 'cox_p_adjusted':
            p_val = gene_info['cox_p_adjusted']
            hr_info = f"HR={gene_info.get('cox_hr_adjusted', 'NA'):.3f}" if pd.notna(gene_info.get('cox_hr_adjusted')) else ""
        elif primary_p_col == 'cox_p_unadjusted':
            p_val = gene_info['cox_p_unadjusted']
            hr_info = f"HR={gene_info.get('cox_hr_unadjusted', 'NA'):.3f}" if pd.notna(gene_info.get('cox_hr_unadjusted')) else ""
        elif primary_p_col == 'logrank_p_value':
            p_val = gene_info['logrank_p_value']
            hr_info = ""
        else:
            p_val = gene_info['p_value']
            hr_info = ""
        
        logger.info(f"  {i+1:2d}. {gene_info['gene']:20s} P={p_val:.2e} {hr_info} "
                   f"(高表达组n={gene_info['high_group_n']}, 低表达组n={gene_info['low_group_n']})")
    
    return results_df

def plot_km_curves_for_top_genes(expr_df, survival_df, results_df, output_dir, top_n=10, logger=None):
    """
    为最显著的基因绘制KM曲线
    
    【分组策略详细说明】
    ===================
    
    1. 基因表达分组方法：
       - 计算每个基因在所有样本中的表达值分布
       - 使用分位数进行分组：
         * 高表达组：表达值 ≥ 60%分位数 (最高40%)
         * 低表达组：表达值 ≤ 40%分位数 (最低40%)
         * 中等表达组：40%分位数 < 表达值 < 60%分位数 (中间20%，排除不参与比较)
    
    2. 分组逻辑原理：
       - 排除中等表达组是为了增强组间对比效果
       - 这种分组方法确保了高低表达组有明显的表达差异
       - 使用40%/60%分位数比传统的25%/75%分位数包含更多样本，提高统计检验功效
       - 避免了任意二分法可能造成的统计偏倚
       
    3. 样本数要求：
       - 每组至少需要5个样本才能进行可靠的生存分析
       - 总样本数至少20个才进行分析
       
    4. 统计检验：
       - Log-rank检验：比较两组生存曲线是否有显著差异
       - 如果有年龄数据，会额外进行Cox回归分析调整年龄因素
       
    5. 结果解释：
       - P值 < 0.05：表示该基因的高低表达对患者生存有显著影响
       - 风险比(HR)：>1表示高表达增加死亡风险，<1表示高表达有保护作用
    
    参数：
        expr_df (DataFrame): 基因表达数据
        survival_df (DataFrame): 生存数据
        results_df (DataFrame): 生存分析结果
        output_dir (Path): 输出目录
        top_n (int): 绘制前N个最显著的基因
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
        
        # 获取基因表达数据
        gene_data = expr_df.loc[gene_name]
        
        # 合并数据
        merged_df = pd.merge(
            pd.DataFrame({'expression': gene_data}),
            survival_df[['time', 'event']],
            left_index=True,
            right_index=True
        )
        
        # 分组（低表达组：最低40%，高表达组：最高40%）
        low_quantile = merged_df['expression'].quantile(0.40)
        high_quantile = merged_df['expression'].quantile(0.60)
        
        high_group = merged_df[merged_df['expression'] >= high_quantile]
        low_group = merged_df[merged_df['expression'] <= low_quantile]
        
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
    
    if logger:
        logger.info(f"✓ 完成所有KM曲线绘制")

def create_summary_plots(results_df, output_dir, logger):
    """
    创建汇总图表
    
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
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    p_values = results_df[p_col].dropna()
    plt.hist(p_values, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(x=0.05, color='red', linestyle='--', label='P=0.05')
    plt.axvline(x=0.01, color='orange', linestyle='--', label='P=0.01')
    plt.xlabel(p_title)
    plt.ylabel('基因数量')
    plt.title(f'{p_title}分布直方图')
    plt.legend()
    
    # 2. -log10(P值) vs 基因排名
    plt.subplot(2, 2, 2)
    ranks = range(1, len(p_values) + 1)
    plt.plot(ranks, -np.log10(p_values), 'o-', markersize=3, alpha=0.7)
    plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='P=0.05')
    plt.axhline(y=-np.log10(0.01), color='orange', linestyle='--', label='P=0.01')
    plt.xlabel('基因排名')
    plt.ylabel(f'-log10({p_title})')
    plt.title('显著性排名图')
    plt.legend()
    
    # 3. 样本数分布
    plt.subplot(2, 2, 3)
    total_samples = results_df['high_group_n'] + results_df['low_group_n']
    plt.hist(total_samples, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    plt.xlabel('分析样本数')
    plt.ylabel('基因数量')
    plt.title('样本数分布')
    
    # 4. 前20个最显著基因的条形图
    plt.subplot(2, 2, 4)
    top_20 = results_df.head(20)
    y_pos = range(len(top_20))
    plt.barh(y_pos, -np.log10(top_20[p_col]), color='coral', alpha=0.8)
    plt.yticks(y_pos, [gene[:15] + '...' if len(gene) > 15 else gene for gene in top_20['gene']])
    plt.xlabel(f'-log10({p_title})')
    plt.title('前20个最显著基因')
    plt.gca().invert_yaxis()
    
    plt.tight_layout()
    
    # 保存汇总图
    summary_path = output_dir / "survival_analysis_summary"
    plt.savefig(f"{summary_path}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{summary_path}.pdf", bbox_inches='tight')
    plt.close()
    
    logger.info(f"✓ 汇总图表已保存到 {summary_path.name}")

def main():
    """
    批量生存分析主函数
    """
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="批量生存分析 - 识别与患者生存预后最相关的基因")
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
        default=20,
        help="绘制前N个最显著基因的KM曲线 (默认: 20)"
    )
    parser.add_argument(
        "--n_cores",
        type=int,
        default=None,
        help="使用的CPU核心数 (默认: CPU总数-1)"
    )
    parser.add_argument(
        "--skip_preview",
        action="store_true",
        help="跳过数据预览确认步骤，直接开始分析"
    )
    parser.add_argument(
        "--adjust_age",
        action="store_true",
        default=True,
        help="是否调整年龄混杂因素 (默认: True)"
    )
    parser.add_argument(
        "--no_adjust_age",
        action="store_true",
        help="不调整年龄混杂因素（使用传统Log-rank检验）"
    )
    
    args = parser.parse_args()
    
    # 创建输出目录
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # 设置日志
    logger = setup_logging(args.output_dir)
    
    # 记录开始时间
    start_time = datetime.now()
    logger.info("="*80)
    logger.info("批量生存分析开始")
    logger.info("="*80)
    logger.info(f"分析开始时间: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"配置参数:")
    logger.info(f"  - 表达数据文件: {args.expression_file.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 临床数据文件: {args.clinical_file.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 输出目录: {args.output_dir.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 绘制KM曲线数: {args.top_n}")
    logger.info(f"  - CPU核心数: {args.n_cores}")
    logger.info(f"  - 跳过数据预览: {'是' if args.skip_preview else '否'}")
    logger.info(f"  - 调整年龄: {'是' if args.adjust_age else '否'}")
    
    try:
        # 1. 准备数据
        expr_df, survival_df = prepare_survival_data(args.expression_file, args.clinical_file, logger)
        if expr_df is None or survival_df is None:
            logger.error("✗ 数据准备失败，终止分析")
            return
        
        # 2. 如果未跳过预览，等待用户确认
        if not args.skip_preview:
            logger.info("="*80)
            logger.info("数据预览文件已保存，请检查以下文件:")
            logger.info("="*80)
            logger.info("1. survival_data_preview.csv - 生存数据预览")
            logger.info("2. expression_data_preview_top10_genes.csv - 表达数据预览")
            logger.info("3. data_summary_statistics.csv - 数据统计信息")
            logger.info("")
            logger.info("请检查这些文件以确认数据质量。")
            logger.info("如果数据看起来正常，请按回车继续分析...")
            logger.info("如果需要停止分析，请按 Ctrl+C")
            
            try:
                input("按回车键继续分析，或按 Ctrl+C 停止...")
                logger.info("✓ 用户确认继续分析")
            except KeyboardInterrupt:
                logger.info("✗ 用户选择停止分析")
                return
        else:
            logger.info("跳过数据预览确认步骤，直接开始分析")
        
        # 3. 批量生存分析
        analysis_start_time = datetime.now()
        logger.info("="*80)
        logger.info(f"正式开始批量生存分析: {analysis_start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("="*80)
        
        results_df = batch_survival_analysis(expr_df, survival_df, logger, args.n_cores)
        if results_df is None:
            logger.error("✗ 批量生存分析失败")
            return
        
        analysis_end_time = datetime.now()
        analysis_duration = (analysis_end_time - analysis_start_time).total_seconds()
        logger.info(f"✓ 分析完成，耗时: {analysis_duration/60:.1f} 分钟")
        
        # 4. 保存结果
        logger.info("="*80)
        logger.info("保存分析结果")
        logger.info("="*80)
        
        results_file = args.output_dir / "all_genes_survival_analysis_results.csv"
        results_df.to_csv(results_file, index=False)
        logger.info(f"✓ 完整结果已保存到: {results_file.name}")
        
        # 5. 保存前100个最显著基因的详细结果
        top_100_file = args.output_dir / "top_100_significant_genes.csv"
        results_df.head(100).to_csv(top_100_file, index=False)
        logger.info(f"✓ 前100个显著基因已保存到: {top_100_file.name}")
        
        # 6. 保存显著基因列表
        # 确定使用哪个P值列
        if 'cox_p_adjusted' in results_df.columns and not results_df['cox_p_adjusted'].isna().all():
            p_col = 'cox_p_adjusted'
        elif 'cox_p_unadjusted' in results_df.columns and not results_df['cox_p_unadjusted'].isna().all():
            p_col = 'cox_p_unadjusted'
        elif 'logrank_p_value' in results_df.columns:
            p_col = 'logrank_p_value'
        else:
            p_col = 'p_value'  # 兼容旧格式
        
        significant_genes = results_df[results_df[p_col] < 0.05]
        if len(significant_genes) > 0:
            sig_file = args.output_dir / "significant_genes_p005.csv"
            significant_genes.to_csv(sig_file, index=False)
            logger.info(f"✓ 显著基因(P<0.05)已保存到: {sig_file.name} (使用{p_col})")
        
        # 7. 绘制KM曲线
        logger.info("="*80)
        logger.info("生成可视化结果")
        logger.info("="*80)
        
        plot_km_curves_for_top_genes(
            expr_df, survival_df, results_df, 
            args.output_dir, args.top_n, logger
        )
        
        # 8. 创建汇总图表
        create_summary_plots(results_df, args.output_dir, logger)
        
        # 记录结束时间
        end_time = datetime.now()
        total_duration = (end_time - start_time).total_seconds()
        
        logger.info("="*80)
        logger.info("批量生存分析完成")
        logger.info("="*80)
        logger.info(f"分析结束时间: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"总耗时: {total_duration/60:.1f} 分钟")
        logger.info(f"平均每个基因耗时: {analysis_duration/len(results_df):.3f} 秒")
        
        # 输出最终结果总结
        logger.info(f"\n最终结果总结:")
        logger.info(f"  - 成功分析基因数: {len(results_df)}")
        logger.info(f"  - 显著基因数(P<0.05): {len(significant_genes)}")
        
        # 计算校正后显著基因数
        if 'cox_p_adjusted_fdr' in results_df.columns:
            adj_sig_count = len(results_df[results_df['cox_p_adjusted_fdr'] < 0.05])
        elif 'cox_p_unadjusted_fdr' in results_df.columns:
            adj_sig_count = len(results_df[results_df['cox_p_unadjusted_fdr'] < 0.05])
        elif 'logrank_p_fdr' in results_df.columns:
            adj_sig_count = len(results_df[results_df['logrank_p_fdr'] < 0.05])
        elif 'p_value_adj' in results_df.columns:
            adj_sig_count = len(results_df[results_df['p_value_adj'] < 0.05])
        else:
            adj_sig_count = 0
        
        logger.info(f"  - 校正后显著基因数: {adj_sig_count}")
        logger.info(f"  - 最小P值: {results_df[p_col].min():.2e}")
        
        # 输出前10个最显著基因
        logger.info(f"\n前10个最显著的预后基因:")
        for i, (_, gene_info) in enumerate(results_df.head(10).iterrows()):
            p_value = gene_info[p_col]
            logger.info(f"  {i+1:2d}. {gene_info['gene']:15s} "
                       f"P={p_value:.2e} "
                       f"(高表达组: {gene_info['high_group_n']}, "
                       f"低表达组: {gene_info['low_group_n']})")
        
        logger.info(f"\n分析成功完成！所有结果已保存到: {args.output_dir.relative_to(PROJECT_ROOT)}")
        
    except Exception as e:
        logger.error(f"✗ 分析过程中发生错误: {e}")
        import traceback
        logger.error(f"错误详情: {traceback.format_exc()}")
        return

if __name__ == "__main__":
    main() 