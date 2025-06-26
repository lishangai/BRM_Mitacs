如果您的数据仅限于 illuminahiseq_rnaseqv2-RSEM_genes_normalized（已标准化的RSEM表达矩阵），以下是 严格基于现有数据 可完成的完整分析步骤（Python实现）：

## 一、数据预处理

加载与基础清洗

import pandas as pd
import numpy as np

加载RSEM标准化数据（行是基因，列是样本）

expr = pd.read_csv("illuminahiseq_rnaseqv2-RSEM_genes_normalized.csv", index_col=0)

检查BRM基因（确保基因符号一致）

print("BRM" in expr.index)  # 若为False，尝试"SMARCA2"（BRM的官方符号）

过滤低表达基因（保留在≥10%样本中TPM>1的基因）

expr = expr[(expr > 1).sum(axis=1) >= 0.1 * expr.shape[1]]

BRM表达分组

取log2(RSEM+1)后分组（中位数）

brm_expr = np.log2(expr.loc["SMARCA2"] + 1)
brm_median = np.median(brm_expr)
groups = np.where(brm_expr > brm_median, "High", "Low")

## 二、差异表达分析（非参数方法）

差异基因筛选

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

deg_results = []
for gene in expr.index:
    high = np.log2(expr.loc[gene, groups == "High"] + 1)
    low = np.log2(expr.loc[gene, groups == "Low"] + 1)
    pval = mannwhitneyu(high, low, alternative="two-sided").pvalue
    logfc = np.mean(high) - np.mean(low)
    deg_results.append([gene, logfc, pval])

deg_df = pd.DataFrame(deg_results, columns=["gene", "log2FC", "pval"])
deg_df["padj"] = fdrcorrection(deg_df["pval"])[1]  # FDR校正
significant_genes = deg_df[deg_df["padj"] < 0.05].sort_values("log2FC", key=abs, ascending=False)

可视化

import seaborn as sns
import matplotlib.pyplot as plt

火山图

plt.figure(figsize=(10,6))
sns.scatterplot(data=deg_df, x="log2FC", y="-log10(padj)", 
                hue=deg_df["padj"] < 0.05, palette={True: "red", False: "gray"})
plt.title("Volcano Plot: BRM-High vs Low")
plt.show()

热图（Top 50 DEGs）

top_genes = significant_genes.head(50)["gene"]
sns.clustermap(np.log2(expr.loc[top_genes] + 1), 
               col_colors=pd.Series(groups, index=expr.columns).map({"High": "red", "Low": "blue"}),
               figsize=(12,8))

## 三、通路分析（仅需基因列表）

GSEA富集

import gseapy

准备排序基因列表（按log2FC）

gene_rank = deg_df.set_index("gene")["log2FC"].sort_values(ascending=False)

运行GSEA（需联网下载基因集）

gsea_results = gseapy.gsea(
    gene_rank=gene_rank,
    gene_sets="KEGG_2021_Human",
    outdir=None
)

可视化前10通路

gseapy.dotplot(gsea_results.results, title="GSEA: Top Pathways", top=10)

Overrepresentation Analysis (ORA)

ora_results = gseapy.enrichr(
    gene_list=significant_genes["gene"].tolist(),
    gene_sets=["GO_Biological_Process_2021"],
    cutoff=0.05
)
ora_results.results.head(10)

四、共表达网络分析
计算BRM相关基因

计算与BRM的Spearman相关系数

corr_results = []
for gene in expr.index:
    rho, pval = spearmanr(expr.loc["SMARCA2"], expr.loc[gene])
    corr_results.append([gene, rho, pval])

corr_df = pd.DataFrame(corr_results, columns=["gene", "rho", "pval"])
corr_df["padj"] = fdrcorrection(corr_df["pval"])[1]
significant_corr = corr_df[(corr_df["padj"] < 0.05) & (abs(corr_df["rho"]) > 0.3)]

可视化

共表达网络（需安装networkx）

import networkx as nx
= nx.Graph()

G.add_node("BRM")
for _, row in significant_corr.iterrows():
    G.add_edge("BRM", row["gene"], weight=row["rho"])

nx.draw(G, with_labels=True, node_size=50, 
        edge_color=["red" if d["weight"]>0 else "blue" for (_,_,d) in G.edges(data=True)])

五、局限性说明
无法完成的分析：

生存分析（无临床数据）

免疫浸润（需反卷积算法支持的原始counts）

甲基化/突变整合（需其他组学数据）
替代建议：

从 cBioPortal 下载匹配的临床数据（样本ID需一致）

使用 ESTIMATE 算法的近似方法（基于特定基因签名）：

          immune_genes = ["CD3D", "CD8A", "GZMB", "IFNG"]
     stromal_genes = ["COL1A1", "COL3A1", "FN1"]
     expr["immune_score"] = expr.loc[immune_genes].mean()
     expr["stromal_score"] = expr.loc[stromal_genes].mean()


六、输出结果清单
分析类型       可交付结果 文件格式

差异表达 DEGs表格（gene, log2FC, padj） CSV
通路富集 GSEA前10通路图表 + ORA结果 PNG + CSV
共表达网络 BRM相关基因网络图 PNG/GraphML
质量控制 分组表达分布箱线图 PNG

