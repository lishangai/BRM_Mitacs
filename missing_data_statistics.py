"""
缺失数据统计脚本
统计CSV文件中每行每列值为0的数据比例
"""

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

def calculate_zero_proportion(data):
    """
    计算数据中值为0的比例
    
    Parameters:
    -----------
    data : pandas.DataFrame
        输入数据
        
    Returns:
    --------
    row_zero_prop : pandas.Series
        每行零值比例
    col_zero_prop : pandas.Series  
        每列零值比例
    """
    # 计算每行零值比例
    row_zero_counts = (data == 0).sum(axis=1)
    row_zero_prop = row_zero_counts / data.shape[1]
    
    # 计算每列零值比例
    col_zero_counts = (data == 0).sum(axis=0)
    col_zero_prop = col_zero_counts / data.shape[0]
    
    return row_zero_prop, col_zero_prop

def plot_zero_proportions(row_zero_prop, col_zero_prop, output_dir="./"):
    """
    绘制零值比例分布图
    
    Parameters:
    -----------
    row_zero_prop : pandas.Series
        每行零值比例
    col_zero_prop : pandas.Series
        每列零值比例
    output_dir : str
        输出目录
    """
    # 创建输出目录
    Path(output_dir).mkdir(exist_ok=True)
    
    # 绘制行零值比例分布
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.hist(row_zero_prop, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('每行零值比例')
    plt.ylabel('频数')
    plt.title('每行零值比例分布')
    plt.grid(True, alpha=0.3)
    
    # 绘制列零值比例分布
    plt.subplot(1, 2, 2)
    plt.hist(col_zero_prop, bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
    plt.xlabel('每列零值比例')
    plt.ylabel('频数')
    plt.title('每列零值比例分布')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/zero_proportion_distribution.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # 绘制热力图（如果数据不太大）
    if len(col_zero_prop) <= 100 and len(row_zero_prop) <= 100:
        plt.figure(figsize=(10, 8))
        
        # 创建一个矩阵来显示前100行和前100列的零值情况
        sample_data = pd.DataFrame({
            'col_zero_prop': col_zero_prop.head(100),
            'row_zero_prop': row_zero_prop.head(100)
        })
        
        plt.subplot(2, 1, 1)
        plt.plot(range(len(col_zero_prop)), col_zero_prop, 'b-', alpha=0.7, linewidth=1)
        plt.xlabel('列索引')
        plt.ylabel('零值比例')
        plt.title('每列零值比例变化趋势')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(2, 1, 2)
        plt.plot(range(len(row_zero_prop)), row_zero_prop, 'r-', alpha=0.7, linewidth=1)
        plt.xlabel('行索引')
        plt.ylabel('零值比例')
        plt.title('每行零值比例变化趋势')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/zero_proportion_trends.png", dpi=300, bbox_inches='tight')
        plt.show()

def generate_summary_report(data, row_zero_prop, col_zero_prop, output_dir="./"):
    """
    生成汇总报告
    
    Parameters:
    -----------
    data : pandas.DataFrame
        原始数据
    row_zero_prop : pandas.Series
        每行零值比例
    col_zero_prop : pandas.Series
        每列零值比例
    output_dir : str
        输出目录
    """
    # 创建输出目录
    Path(output_dir).mkdir(exist_ok=True)
    
    # 生成汇总统计
    summary = {
        "数据基本信息": {
            "总行数": data.shape[0],
            "总列数": data.shape[1],
            "总元素数": data.size,
            "零值总数": (data == 0).sum().sum(),
            "整体零值比例": (data == 0).sum().sum() / data.size
        },
        "行零值比例统计": {
            "最小值": row_zero_prop.min(),
            "最大值": row_zero_prop.max(),
            "平均值": row_zero_prop.mean(),
            "中位数": row_zero_prop.median(),
            "标准差": row_zero_prop.std(),
            "完全为零的行数": (row_zero_prop == 1.0).sum(),
            "完全非零的行数": (row_zero_prop == 0.0).sum()
        },
        "列零值比例统计": {
            "最小值": col_zero_prop.min(),
            "最大值": col_zero_prop.max(),
            "平均值": col_zero_prop.mean(),
            "中位数": col_zero_prop.median(),
            "标准差": col_zero_prop.std(),
            "完全为零的列数": (col_zero_prop == 1.0).sum(),
            "完全非零的列数": (col_zero_prop == 0.0).sum()
        }
    }
    
    # 打印汇总报告
    print("\n" + "="*60)
    print("           缺失数据（零值）统计报告")
    print("="*60)
    
    for category, stats in summary.items():
        print(f"\n{category}:")
        print("-" * 40)
        for key, value in stats.items():
            if isinstance(value, float):
                print(f"  {key}: {value:.4f}")
            else:
                print(f"  {key}: {value}")
    
    # 保存详细结果到CSV
    row_results = pd.DataFrame({
        'row_index': range(len(row_zero_prop)),
        'zero_proportion': row_zero_prop,
        'zero_count': (data == 0).sum(axis=1),
        'total_count': data.shape[1]
    })
    row_results.to_csv(f"{output_dir}/row_zero_statistics.csv", index=False, encoding='utf-8-sig')
    
    col_results = pd.DataFrame({
        'column_name': data.columns if hasattr(data, 'columns') else range(len(col_zero_prop)),
        'zero_proportion': col_zero_prop,
        'zero_count': (data == 0).sum(axis=0),
        'total_count': data.shape[0]
    })
    col_results.to_csv(f"{output_dir}/column_zero_statistics.csv", index=False, encoding='utf-8-sig')
    
    # 保存汇总报告
    with open(f"{output_dir}/summary_report.txt", 'w', encoding='utf-8') as f:
        f.write("缺失数据（零值）统计报告\n")
        f.write("="*60 + "\n\n")
        
        for category, stats in summary.items():
            f.write(f"{category}:\n")
            f.write("-" * 40 + "\n")
            for key, value in stats.items():
                if isinstance(value, float):
                    f.write(f"  {key}: {value:.4f}\n")
                else:
                    f.write(f"  {key}: {value}\n")
            f.write("\n")
    
    print(f"\n详细结果已保存到：")
    print(f"  - 行统计: {output_dir}/row_zero_statistics.csv")
    print(f"  - 列统计: {output_dir}/column_zero_statistics.csv") 
    print(f"  - 汇总报告: {output_dir}/summary_report.txt")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='统计CSV文件中每行每列零值比例')
    parser.add_argument('input_file', help='输入CSV文件路径')
    parser.add_argument('--output_dir', '-o', default='./', help='输出目录（默认当前目录）')
    parser.add_argument('--delimiter', '-d', default=',', help='CSV分隔符（默认逗号）')
    parser.add_argument('--header', action='store_true', help='CSV文件是否包含表头')
    parser.add_argument('--encoding', '-e', default='utf-8', help='文件编码（默认utf-8）')
    parser.add_argument('--sample_size', '-s', type=int, help='样本大小（用于大文件采样）')
    
    args = parser.parse_args()
    
    try:
        print(f"正在读取文件: {args.input_file}")
        
        # 读取CSV文件
        read_kwargs = {
            'delimiter': args.delimiter,
            'encoding': args.encoding,
            'header': 0 if args.header else None
        }
        
        if args.sample_size:
            # 对大文件进行采样
            data = pd.read_csv(args.input_file, nrows=args.sample_size, **read_kwargs)
            print(f"采样读取前 {args.sample_size} 行数据")
        else:
            data = pd.read_csv(args.input_file, **read_kwargs)
        
        print(f"数据形状: {data.shape}")
        
        # 确保数据是数值型
        numeric_columns = data.select_dtypes(include=[np.number]).columns
        if len(numeric_columns) < data.shape[1]:
            print(f"注意: 检测到非数值列，将只分析数值列 ({len(numeric_columns)}/{data.shape[1]})")
            data = data[numeric_columns]
        
        # 计算零值比例
        print("正在计算零值比例...")
        row_zero_prop, col_zero_prop = calculate_zero_proportion(data)
        
        # 生成汇总报告
        generate_summary_report(data, row_zero_prop, col_zero_prop, args.output_dir)
        
        # 绘制图表
        print("正在生成可视化图表...")
        plot_zero_proportions(row_zero_prop, col_zero_prop, args.output_dir)
        
        print("\n分析完成！")
        
    except Exception as e:
        print(f"错误: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main()) 