"""
GO富集分析脚本
==============

功能概述：
    这个脚本用于对显著差异表达基因进行GO（Gene Ontology）富集分析，
    识别在生物学过程、细胞组分和分子功能方面富集的基因集合。

主要功能：
    1. 读取显著差异表达基因列表
    2. 使用Enrichr数据库进行GO富集分析
    3. 生成富集分析结果的统计表格
    4. 创建富集分析结果的可视化图表

输入：
    - 显著差异表达基因文件（CSV格式）
    
输出：
    - GO富集分析结果表格（CSV文件）
    - 富集分析点图可视化（PNG和PDF格式）
    - 详细的分析日志（log.txt）

作者：BRM_Mitacs项目组
版本：1.2
更新：添加基因ID映射功能，处理?|ID格式的基因条目
"""

import argparse
from pathlib import Path
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import logging
import sys
from datetime import datetime
import traceback
import re # Added for gene recognition statistics

# 定义项目根目录路径
PROJECT_ROOT = Path(__file__).resolve().parents[1]

# 基因ID到基因名称的映射字典
# 用于处理以?|开头的基因条目，将Entrez Gene ID映射到正确的基因符号
GENE_ID_MAPPING = {
    "388795": "EFCAB8",      # EFCAB8 基因
    "57714": "KIAA1618",     # KIAA1618 基因  
    "10431": "TIMM23",       # TIMM23 基因
    # 可以根据需要添加更多映射关系
    # "基因ID": "基因符号",
}

def map_gene_symbols(gene_list, logger):
    """
    映射基因符号，处理?|ID格式的基因条目
    
    参数：
        gene_list (list): 原始基因列表
        logger: 日志记录器
        
    返回：
        list: 映射后的基因列表
    """
    logger.info("开始基因符号映射处理...")
    
    mapped_genes = []
    unmapped_count = 0
    mapped_count = 0
    
    for gene in gene_list:
        if gene.startswith("?|"):
            # 提取基因ID
            gene_id = gene.split("|")[1]
            
            # 查找映射
            if gene_id in GENE_ID_MAPPING:
                mapped_gene = GENE_ID_MAPPING[gene_id]
                mapped_genes.append(mapped_gene)
                mapped_count += 1
                logger.info(f"  映射: {gene} → {mapped_gene}")
            else:
                # 如果没有找到映射，保留原始格式但记录警告
                mapped_genes.append(gene)
                unmapped_count += 1
                logger.warning(f"  未找到映射: {gene} (ID: {gene_id})")
        else:
            # 正常的基因符号，直接添加
            mapped_genes.append(gene)
    
    logger.info(f"基因符号映射统计:")
    logger.info(f"  - 总基因数: {len(gene_list)}")
    logger.info(f"  - 成功映射: {mapped_count}")
    logger.info(f"  - 未映射的?|ID格式: {unmapped_count}")
    logger.info(f"  - 最终基因数: {len(mapped_genes)}")
    
    if unmapped_count > 0:
        logger.warning(f"发现 {unmapped_count} 个未映射的基因ID。")
        logger.warning("您可以在脚本顶部的 GENE_ID_MAPPING 字典中添加更多映射关系。")
        logger.warning("格式: \"基因ID\": \"基因符号\"")
    
    return mapped_genes

def setup_logging(output_dir):
    """
    设置日志记录系统
    
    参数：
        output_dir (Path): 日志文件保存目录
        
    返回：
        logger: 配置好的日志记录器
    """
    # 创建日志文件路径
    log_file = output_dir / "log.txt"
    
    # 创建日志记录器
    logger = logging.getLogger('GO_Enrichment')
    logger.setLevel(logging.INFO)
    
    # 清除现有的处理器（避免重复记录）
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 创建文件处理器
    file_handler = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.INFO)
    
    # 创建控制台处理器
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    
    # 创建格式化器
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # 添加处理器到日志记录器
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

def log_system_info(logger):
    """
    记录系统和环境信息
    
    参数：
        logger: 日志记录器
    """
    logger.info("=" * 80)
    logger.info("GO富集分析系统信息")
    logger.info("=" * 80)
    logger.info(f"分析开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Python版本: {sys.version}")
    logger.info(f"工作目录: {PROJECT_ROOT}")
    
    # 尝试获取关键库的版本信息
    try:
        import pandas as pd
        logger.info(f"Pandas版本: {pd.__version__}")
    except:
        pass
    
    try:
        import matplotlib
        logger.info(f"Matplotlib版本: {matplotlib.__version__}")
    except:
        pass
    
    try:
        import gseapy
        logger.info(f"GSEApy版本: {gseapy.__version__}")
    except:
        pass
    
    logger.info("=" * 80)

def run_go_enrichment(gene_list, output_dir, analysis_name, logger):
    """
    执行GO富集分析的核心函数
    
    功能说明：
        使用gseapy库调用Enrichr在线数据库，对输入的基因列表进行GO富集分析。
        分析三个主要的GO数据库：生物学过程、细胞组分、分子功能。
        
    参数：
        gene_list (list): 基因符号列表，用于富集分析
        output_dir (Path): 结果保存目录的路径对象
        analysis_name (str): 分析名称，用于创建子目录
        logger: 日志记录器
        
    输出：
        - 每个GO数据库的富集结果CSV文件
        - 每个GO数据库前10个最显著条目的点图可视化
    """
    logger.info("开始GO富集分析...")
    logger.info("-" * 60)
    
    # 定义要分析的GO基因集合数据库
    # 使用2021年版本的GO数据库，包含最新的基因注释信息
    gene_sets = [
        'GO_Biological_Process_2021',    # GO生物学过程数据库
        'GO_Cellular_Component_2021',    # GO细胞组分数据库
        'GO_Molecular_Function_2021'     # GO分子功能数据库
    ]

    logger.info(f"分析配置:")
    logger.info(f"  - 将要分析的基因集合: {', '.join(gene_sets)}")
    logger.info(f"  - 用于分析的基因数量: {len(gene_list)}")
    logger.info(f"  - P值阈值: 0.05")
    logger.info(f"  - 物种: 人类 (Homo sapiens)")
    
    # 记录输入基因列表的详细信息
    logger.info(f"\n输入基因列表详情:")
    logger.info(f"  - 基因总数: {len(gene_list)}")
    if len(gene_list) <= 20:
        logger.info(f"  - 基因列表: {', '.join(gene_list)}")
    else:
        logger.info(f"  - 前10个基因: {', '.join(gene_list[:10])}")
        logger.info(f"  - 后10个基因: {', '.join(gene_list[-10:])}")

    # 使用Enrichr进行富集分析
    # Enrichr是一个流行的在线基因富集分析工具
    logger.info(f"\n正在连接Enrichr服务器...")
    try:
        start_time = datetime.now()
        logger.info(f"  - 开始时间: {start_time.strftime('%H:%M:%S')}")
        
        enr = gp.enrichr(
            gene_list=gene_list,                              # 输入基因列表
            gene_sets=gene_sets,                              # GO数据库列表
            organism='human',                                 # 物种：人类
            outdir=str(output_dir / "enrichr_temp"),         # 临时文件目录
            cutoff=0.05                                       # P值阈值，只保留显著结果
        )
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        logger.info(f"  - 完成时间: {end_time.strftime('%H:%M:%S')}")
        logger.info(f"  - 分析耗时: {duration:.2f} 秒")
        logger.info(f"✓ 富集分析成功完成")
        
        # 统计：有多少输入基因被数据库识别（出现在任一条目的 Genes 列中）
        try:
            recognized_genes = set()
            if hasattr(enr, 'results') and enr.results is not None and 'Genes' in enr.results.columns:
                for gstr in enr.results['Genes'].dropna().astype(str):
                    # Enrichr 的 Genes 字段通常以 ';' 分隔，也可能为 ','
                    parts = [p.strip().upper() for p in re.split(r'[;,]', gstr) if p and p.strip()]
                    recognized_genes.update(parts)
            input_genes_set = set([g.upper() for g in gene_list])
            recognized_from_input = input_genes_set.intersection(recognized_genes)
            not_recognized = sorted(list(input_genes_set - recognized_from_input))
            logger.info("输入基因识别统计（基于Enrichr返回的Genes字段）:")
            logger.info(f"  - 输入基因数: {len(input_genes_set)}")
            logger.info(f"  - 被数据库成功识别: {len(recognized_from_input)}")
            logger.info(f"  - 未被识别: {len(not_recognized)}")
            if len(not_recognized) <= 20 and len(not_recognized) > 0:
                logger.info(f"  - 未识别列表: {', '.join(not_recognized)}")
            elif len(not_recognized) > 20:
                logger.info(f"  - 未识别前10例: {', '.join(not_recognized[:10])}")
        except Exception as e2:
            logger.warning(f"无法统计数据库识别基因数: {e2}")
        
    except Exception as e:
        logger.error(f"✗ 富集分析过程中发生错误: {e}")
        logger.error(f"错误详情: {traceback.format_exc()}")
        logger.error("请检查以下问题:")
        logger.error("  1. 网络连接是否正常")
        logger.error("  2. Enrichr服务器是否可用")
        logger.error("  3. 输入基因格式是否正确")
        return

    logger.info(f"\n开始处理富集分析结果...")

    # 为当前分析创建专门的结果目录
    results_path = output_dir / analysis_name
    results_path.mkdir(parents=True, exist_ok=True)
    logger.info(f"结果保存目录: {results_path.relative_to(PROJECT_ROOT)}")

    # 总体统计信息
    total_results = len(enr.results) if hasattr(enr, 'results') and enr.results is not None else 0
    logger.info(f"\n总体富集分析统计:")
    logger.info(f"  - 总结果条目数: {total_results}")

    # 处理并保存每个基因集合的结果
    analysis_summary = {
        'gene_sets': {},
        'total_significant_terms': 0,
        'files_generated': []
    }
    
    for i, gene_set_name in enumerate(gene_sets):
        logger.info(f"\n{'='*50}")
        logger.info(f"正在处理: {gene_set_name}")
        logger.info(f"{'='*50}")
        
        # 从完整结果中提取当前基因集合的数据
        result_df = enr.results[enr.results['Gene_set'] == gene_set_name]
        
        logger.info(f"原始结果条目数: {len(result_df)}")

        # 检查是否有显著的富集结果
        if result_df.empty:
            logger.warning(f"⚠ {gene_set_name} 中未发现显著富集的条目")
            analysis_summary['gene_sets'][gene_set_name] = {
                'significant_terms': 0,
                'files_created': []
            }
            continue

        # 详细统计信息
        significant_count = len(result_df)
        analysis_summary['total_significant_terms'] += significant_count
        analysis_summary['gene_sets'][gene_set_name] = {
            'significant_terms': significant_count,
            'files_created': []
        }
        
        logger.info(f"✓ 发现 {significant_count} 个显著富集条目")
        
        # 统计P值分布
        if 'P-value' in result_df.columns:
            p_values = result_df['P-value']
            logger.info(f"P值统计:")
            logger.info(f"  - 最小P值: {p_values.min():.2e}")
            logger.info(f"  - 最大P值: {p_values.max():.2e}")
            logger.info(f"  - 平均P值: {p_values.mean():.2e}")
            logger.info(f"  - P值 < 0.001 的条目数: {(p_values < 0.001).sum()}")
            logger.info(f"  - P值 < 0.01 的条目数: {(p_values < 0.01).sum()}")

        # 保存完整的富集分析结果到CSV文件
        csv_path = results_path / f"{gene_set_name}_results.csv"
        try:
            result_df.to_csv(csv_path, index=False)
            file_size = csv_path.stat().st_size / 1024  # KB
            logger.info(f"✓ 已保存完整结果到: {csv_path.name} ({file_size:.1f} KB)")
            analysis_summary['gene_sets'][gene_set_name]['files_created'].append(csv_path.name)
            analysis_summary['files_generated'].append(csv_path.name)
            
            # 记录列信息
            logger.info(f"结果表格包含的列: {', '.join(result_df.columns)}")
            
        except Exception as e:
            logger.error(f"✗ 保存CSV文件时发生错误: {e}")

        # 为最显著的前10个条目创建点图可视化
        # 按P值排序以获得最显著的结果
        top_terms = result_df.sort_values(by='P-value').head(10)
        
        logger.info(f"\n准备生成可视化图表:")
        logger.info(f"  - 用于绘图的条目数: {len(top_terms)}")
        
        if len(top_terms) > 0:
            # 显示最显著的前5个条目
            logger.info(f"前5个最显著的富集条目:")
            for idx, (_, row) in enumerate(top_terms.head(5).iterrows(), 1):
                term_name = row.get('Term', 'Unknown')[:50]  # 限制长度
                p_value = row.get('P-value', 'N/A')
                logger.info(f"  {idx}. {term_name}... (P={p_value:.2e})")
        
        if not top_terms.empty:
            try:
                logger.info(f"正在生成 {gene_set_name} 的点图可视化...")
                
                # 使用gseapy内置的点图绘制功能
                fig, ax = plt.subplots(figsize=(12, 8))
                gp.plot.dotplot(top_terms, title=f"{gene_set_name} - Top 10 Enriched Terms", ax=ax, ofname=None)
                
                # 保存图表为PNG和PDF格式
                plot_path = results_path / f"{gene_set_name}_dotplot"
                
                # 保存PNG格式
                png_path = f"{plot_path}.png"
                fig.savefig(png_path, dpi=300, bbox_inches='tight')
                png_size = Path(png_path).stat().st_size / 1024  # KB
                logger.info(f"✓ 保存PNG图表: {plot_path.name}.png ({png_size:.1f} KB)")
                
                # 保存PDF格式
                pdf_path = f"{plot_path}.pdf"
                fig.savefig(pdf_path, bbox_inches='tight')
                pdf_size = Path(pdf_path).stat().st_size / 1024  # KB
                logger.info(f"✓ 保存PDF图表: {plot_path.name}.pdf ({pdf_size:.1f} KB)")
                
                plt.close(fig)  # 关闭图形以释放内存
                
                analysis_summary['gene_sets'][gene_set_name]['files_created'].extend([
                    f"{plot_path.name}.png", 
                    f"{plot_path.name}.pdf"
                ])
                analysis_summary['files_generated'].extend([
                    f"{plot_path.name}.png", 
                    f"{plot_path.name}.pdf"
                ])

            except Exception as e:
                logger.error(f"✗ 无法为 {gene_set_name} 生成图表: {e}")
                logger.error(f"错误详情: {traceback.format_exc()}")
        else:
            logger.warning(f"⚠ {gene_set_name} 没有足够的条目用于绘图")

    # 输出分析总结
    logger.info(f"\n" + "="*80)
    logger.info("GO富集分析完成总结")
    logger.info("="*80)
    logger.info(f"总体统计:")
    logger.info(f"  - 分析的基因数量: {len(gene_list)}")
    logger.info(f"  - 分析的GO数据库数量: {len(gene_sets)}")
    logger.info(f"  - 总显著富集条目数: {analysis_summary['total_significant_terms']}")
    logger.info(f"  - 生成的文件总数: {len(analysis_summary['files_generated'])}")
    
    logger.info(f"\n各数据库详细结果:")
    for gene_set, stats in analysis_summary['gene_sets'].items():
        logger.info(f"  {gene_set}:")
        logger.info(f"    - 显著条目数: {stats['significant_terms']}")
        logger.info(f"    - 生成文件数: {len(stats['files_created'])}")
        if stats['files_created']:
            logger.info(f"    - 文件列表: {', '.join(stats['files_created'])}")
    
    return analysis_summary


def main():
    """
    主函数：GO富集分析工作流程
    
    工作流程：
        1. 解析命令行参数
        2. 设置日志系统
        3. 创建输出目录
        4. 加载显著差异表达基因列表
        5. 数据预处理（基因符号标准化）
        6. 执行GO富集分析
        7. 生成结果报告
    """
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="对显著差异表达基因列表执行GO富集分析")
    
    # 输入文件参数：显著差异表达基因文件路径
    parser.add_argument(
        "--input_file",
        type=Path,
        default=PROJECT_ROOT / "results" / "brm_analysis_csv" / "significant_DEGs.csv",
        help="显著差异表达基因文件的路径。默认: results/brm_analysis_csv/significant_DEGs.csv"
    )
    
    # 输出目录参数：分析结果保存目录
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=PROJECT_ROOT / "results" / "go_enrichment",
        help="分析结果保存目录。默认: results/go_enrichment"
    )
    
    # 分析名称参数：用于创建输出子目录
    parser.add_argument(
        "--analysis_name",
        type=str,
        default="SMARCA2_analysis",
        help="当前分析运行的名称，用于输出子目录的命名。"
    )
    
    # 解析命令行参数
    args = parser.parse_args()

    # 创建输出目录（如果不存在）
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # 设置日志系统
    logger = setup_logging(args.output_dir)
    
    # 记录系统信息
    log_system_info(logger)

    # 打印分析配置信息
    logger.info("GO富集分析工作流程配置")
    logger.info("-" * 60)
    logger.info(f"配置参数:")
    logger.info(f"  - 输入文件      : {args.input_file.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 输出目录      : {args.output_dir.relative_to(PROJECT_ROOT)}")
    logger.info(f"  - 分析名称      : {args.analysis_name}")
    logger.info(f"  - 日志文件      : {(args.output_dir / 'log.txt').relative_to(PROJECT_ROOT)}")

    # 加载和预处理基因列表
    try:
        logger.info(f"\n开始加载显著差异表达基因列表...")
        logger.info(f"输入文件路径: {args.input_file}")
        
        # 检查文件是否存在
        if not args.input_file.exists():
            raise FileNotFoundError(f"输入文件不存在: {args.input_file}")
        
        # 获取文件大小
        file_size = args.input_file.stat().st_size / 1024  # KB
        logger.info(f"文件大小: {file_size:.1f} KB")
        
        # 读取CSV文件
        logger.info("正在读取CSV文件...")
        deg_df = pd.read_csv(args.input_file)
        logger.info(f"✓ 成功读取文件，包含 {len(deg_df)} 行数据")
        logger.info(f"文件列名: {', '.join(deg_df.columns.tolist())}")
        
        # 验证文件格式：确保包含'gene'列
        if 'gene' not in deg_df.columns:
            raise ValueError("输入文件必须包含'gene'列。")
        
        logger.info(f"开始数据预处理...")
        logger.info(f"原始基因条目数: {len(deg_df)}")
        
        # 数据预处理：
        # 1. 分割基因符号（处理类似'GENE|ID'格式的基因名）
        # 2. 转换为大写字母（标准化基因符号格式）
        # 3. 去除缺失值
        # 4. 去重，获取唯一基因列表
        
        # 检查基因符号格式
        sample_genes = deg_df['gene'].head(10).tolist()
        logger.info(f"样本基因符号: {sample_genes}")
        
        # 智能处理基因符号：先映射?|ID格式，再处理其他格式
        raw_gene_list = deg_df['gene'].tolist()
        logger.info(f"原始基因条目数: {len(raw_gene_list)}")
        
        # 第一步：应用基因ID映射（处理?|ID格式）
        mapped_gene_list = map_gene_symbols(raw_gene_list, logger)
        
        # 第二步：处理正常的基因符号（分割|并转大写）
        processed_genes = []
        for gene in mapped_gene_list:
            if '|' in gene and not gene.startswith('?|'):
                # 对于正常的 GENE|ID 格式，取基因符号部分
                gene_symbol = gene.split('|')[0].upper()
                processed_genes.append(gene_symbol)
            else:
                # 对于已映射的基因或单独的基因符号，直接使用
                processed_genes.append(gene.upper())
        
        logger.info(f"处理后基因数: {len(processed_genes)}")
        
        # 去除缺失值和空字符串
        no_na_genes = [gene for gene in processed_genes if gene and gene.strip() and gene != '?']
        logger.info(f"去除缺失值和无效符号后: {len(no_na_genes)}")
        
        # 去重
        gene_list = list(set(no_na_genes))
        logger.info(f"去重后最终基因数: {len(gene_list)}")
        
        # 统计信息
        logger.info(f"\n基因预处理统计:")
        logger.info(f"  - 原始条目数: {len(deg_df)}")
        logger.info(f"  - 处理后有效基因数: {len(no_na_genes)}")
        logger.info(f"  - 重复基因数: {len(no_na_genes) - len(gene_list)}")
        logger.info(f"  - 最终唯一基因数: {len(gene_list)}")
        
        # 显示基因列表样本
        if len(gene_list) <= 10:
            logger.info(f"最终基因列表: {', '.join(gene_list)}")
        else:
            logger.info(f"基因列表样本 (前10个): {', '.join(gene_list[:10])}")
        
    except FileNotFoundError:
        logger.error(f"✗ 错误: 在路径 {args.input_file} 未找到输入文件")
        return
    except ValueError as e:
        logger.error(f"✗ 处理输入文件时发生错误: {e}")
        return
    except Exception as e:
        logger.error(f"✗ 读取文件时发生未知错误: {e}")
        logger.error(f"错误详情: {traceback.format_exc()}")
        return

    # 验证基因列表不为空
    if not gene_list:
        logger.error("✗ 输入文件中未找到有效基因。终止分析。")
        return

    # 执行GO富集分析
    logger.info(f"\n" + "="*80)
    logger.info("开始执行GO富集分析")
    logger.info("="*80)
    
    try:
        analysis_start_time = datetime.now()
        logger.info(f"分析开始时间: {analysis_start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        summary = run_go_enrichment(gene_list, args.output_dir, args.analysis_name, logger)
        
        analysis_end_time = datetime.now()
        total_duration = (analysis_end_time - analysis_start_time).total_seconds()
        
        logger.info(f"\n" + "="*80)
        logger.info("GO富集分析流程完成")
        logger.info("="*80)
        logger.info(f"分析结束时间: {analysis_end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"总耗时: {total_duration:.2f} 秒 ({total_duration/60:.1f} 分钟)")
        logger.info(f"所有结果已保存到: {args.output_dir.relative_to(PROJECT_ROOT) / args.analysis_name}")
        
        logger.info(f"\n生成的文件类型:")
        logger.info(f"  - CSV结果文件: 包含完整的富集分析结果表格")
        logger.info(f"  - PNG图表文件: 高分辨率(300 DPI)的富集结果可视化")
        logger.info(f"  - PDF图表文件: 矢量格式的富集结果可视化")
        logger.info(f"  - 日志文件: 详细的分析过程记录")
        
        logger.info(f"\n分析成功完成！")
        
    except Exception as e:
        logger.error(f"✗ GO富集分析过程中发生错误: {e}")
        logger.error(f"错误详情: {traceback.format_exc()}")
        return

# 脚本入口点
if __name__ == "__main__":
    main() 