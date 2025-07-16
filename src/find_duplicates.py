import pandas as pd
from pathlib import Path

# 读取文件
PROJECT_ROOT = Path.cwd()
input_file = PROJECT_ROOT / "results" / "brm_analysis_csv" / "significant_DEGs.csv"

print("正在分析重复基因...")
print(f"文件路径: {input_file}")

# 读取数据
deg_df = pd.read_csv(input_file)
print(f"总基因条目数: {len(deg_df)}")

# 获取原始基因列表
original_genes = deg_df['gene'].tolist()
print(f"原始基因列表长度: {len(original_genes)}")

# 处理基因符号（与go_enrichment.py中相同的处理方式）
processed_genes = deg_df['gene'].str.split('|').str[0].str.upper()
print(f"处理后基因列表长度: {len(processed_genes)}")

# 去除缺失值
no_na_genes = processed_genes.dropna()
print(f"去除缺失值后长度: {len(no_na_genes)}")

# 找出重复的基因
print("\n分析重复基因...")

# 创建基因计数字典
gene_counts = {}
gene_indices = {}

for i, gene in enumerate(no_na_genes):
    if gene in gene_counts:
        gene_counts[gene] += 1
        gene_indices[gene].append(i)
    else:
        gene_counts[gene] = 1
        gene_indices[gene] = [i]

# 找出重复的基因
duplicated_genes = {gene: count for gene, count in gene_counts.items() if count > 1}

print(f"发现 {len(duplicated_genes)} 个重复的基因:")
print(f"重复基因总条目数: {sum(duplicated_genes.values())}")
print(f"去重后减少的条目数: {sum(duplicated_genes.values()) - len(duplicated_genes)}")

print("\n重复基因详细信息:")
print("="*80)

for gene, count in duplicated_genes.items():
    print(f"\n重复基因: {gene} (出现 {count} 次)")
    indices = gene_indices[gene]
    
    for i, idx in enumerate(indices):
        original_name = original_genes[idx]
        log2fc = deg_df.iloc[idx]['log2FC']
        pval = deg_df.iloc[idx]['pval']
        padj = deg_df.iloc[idx]['padj']
        
        print(f"  第 {i+1} 次出现:")
        print(f"    - 行号: {idx + 2}")  # +2 because of header and 0-based indexing
        print(f"    - 原始名称: {original_name}")
        print(f"    - log2FC: {log2fc:.4f}")
        print(f"    - P值: {pval:.2e}")
        print(f"    - 调整后P值: {padj:.2e}")

print("\n" + "="*80)
print("总结:")
print(f"- 原始基因条目: {len(deg_df)}")
print(f"- 处理后基因条目: {len(no_na_genes)}")
print(f"- 唯一基因数: {len(no_na_genes) - sum(duplicated_genes.values()) + len(duplicated_genes)}")
print(f"- 重复基因种类: {len(duplicated_genes)}")
print(f"- 去除的重复条目: {sum(duplicated_genes.values()) - len(duplicated_genes)}") 