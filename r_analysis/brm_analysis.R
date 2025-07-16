## ---- 1. 环境设置与数据加载 ----

# 设置CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# # 强制重装可能损坏的plyr包 (已在首次运行时修复，现已注释掉)
# install.packages("plyr")
# # 强制重装可能损坏的reshape2包 (已在首次运行时修复，现已注释掉)
# install.packages("reshape2")

# # 强制安装/重装 org.Hs.eg.db (已在首次运行时修复，现已注释掉)
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# BiocManager::install("org.Hs.eg.db")

# Function to check and install packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      BiocManager::install(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# 检查 BiocManager 是否已安装
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 定义需要的包
required_packages <- c(
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "clusterProfiler",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "enrichplot",
  "dplyr",
  "tibble",
  "ggrepel",
  "viridis"
)

# 安装/加载包
install_if_missing(required_packages)

# 确保核心注释库已安装
if (!require("org.Hs.eg.db", character.only = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
  library("org.Hs.eg.db", character.only = TRUE)
}

# 定义文件路径
counts_file <- "OV.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"
# phenotype_file <- "../OV.clin.merged.txt" # 已禁用，不再需要外部临床文件

# 创建输出目录
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

cat("开始加载和解析TCGA Level 3 RSEM文件...\n")

# ---- 2. TCGA RSEM 数据解析 ----

# 读取第一行获取样本ID
header <- read.delim(counts_file, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
# 从第2列开始，每隔3列是一个新的样本ID (raw_count, scaled_estimate, transcript_id)
sample_ids_raw <- as.character(header[1, seq(2, ncol(header), by = 3)])

# 读取数据主体，跳过第一行
raw_data <- read.delim(counts_file, skip = 1, header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# 提取我们需要的 raw_count 列
# raw_count 列是第1, 4, 7, ... 列
raw_count_columns <- seq(1, ncol(raw_data), by = 3)
counts_data <- raw_data[, raw_count_columns]

# 将列名设置为原始样本ID
colnames(counts_data) <- sample_ids_raw

# -- 正确的基因ID处理与去重流程 (保留原始ID) --

# 1. 生成原始基因ID向量 (不做任何改动)
gene_ids_raw <- rownames(raw_data)

# 2. 将原始基因ID作为一列添加到数据中，用于后续聚合
counts_data$gene_id_raw <- gene_ids_raw

# 3. 移除任何以 '?' 开头的基因
counts_data <- counts_data[!grepl("^\\?", counts_data$gene_id_raw), ]

# 4. 使用aggregate函数按原始基因ID分组，并计算均值
cat("正在合并重复基因的表达量...\n")
counts_data_agg <- aggregate(. ~ gene_id_raw, data = counts_data, FUN = mean)

# 5. 将聚合后的、唯一的原始基因ID设置为行名
rownames(counts_data_agg) <- counts_data_agg$gene_id_raw

# 6. 移除临时的gene_id_raw列
counts_data_agg$gene_id_raw <- NULL
counts_data <- counts_data_agg

# 确保所有计数都是整数
counts_data <- round(counts_data)

cat(paste("成功解析数据，得到", nrow(counts_data), "个基因和", ncol(counts_data), "个样本。\n"))


# ---- 3. 动态创建最小化样本信息表 (Phenotype) ----

cat("正在动态创建最小化样本信息表以满足DESeq2要求...\n")
    
# 直接从counts矩阵的列名（即样本ID）创建样本信息表
phenotype_data <- data.frame(row.names = colnames(counts_data))
    
cat(paste("成功创建样本信息表，包含", nrow(phenotype_data), "个样本。\n"))
# 数据已经天然对齐，无需进行复杂的样本匹配。


# ---- 4. DESeq2 分析 ----

# 提取 BRM (SMARCA2) 的表达量 (需要用包含Entrez ID的原始基因名)
brm_gene_id_raw <- rownames(counts_data)[grepl("^SMARCA2\\|", rownames(counts_data))]

if(length(brm_gene_id_raw) == 1) {
    # 将SMARCA2表达量添加到我们动态创建的phenotype_data中
    # 根据DESeq2的建议，我们对这个连续变量进行中心化和标准化(scaling)
    phenotype_data$SMARCA2_expression <- scale(as.numeric(counts_data[brm_gene_id_raw, rownames(phenotype_data)]))
    # 从counts_data中移除BRM基因
    counts_data <- counts_data[rownames(counts_data) != brm_gene_id_raw, ]
} else {
    stop("BRM (SMARCA2) gene ID not found or found multiple matches in the counts data.")
}

# 过滤低表达基因
keep <- rowSums(counts_data >= 10) >= (ncol(counts_data) * 0.1)
counts_data <- counts_data[keep,]

cat(paste("过滤低表达基因后，保留", nrow(counts_data), "个基因进行分析。\n"))

# 创建DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = phenotype_data,
                              design = ~ SMARCA2_expression)

# ---- 4.1. 质量控制 (QC) 与PCA分析 ----
cat("正在进行质量控制(QC)与PCA分析...\n")

# 使用vst进行方差稳定变换，适用于中等到大型数据集
vsd <- vst(dds, blind = FALSE)

# 绘制PCA图，并根据SMARCA2表达量进行着色
pca_plot <- plotPCA(vsd, intgroup = "SMARCA2_expression") +
  theme_bw(base_size = 16) + # 使用更清晰的学术主题
  labs(
    title = "样本主成分分析 (PCA)",
    subtitle = "按SMARCA2表达量着色",
    x = "PC1",
    y = "PC2",
    color = "SMARCA2 Expression\n(Scaled)" # 修改图例标题
  ) +
  scale_color_viridis_c() + # 使用美观且色盲友好的viridis连续色谱
  geom_point(size = 4, alpha = 0.8) + # 让点更突出
  theme(legend.position = "right")

ggsave("plots/PCA_by_SMARCA2_expression.png", plot = pca_plot, width = 10, height = 8, dpi = 300)
cat("PCA图已保存到 r_analysis/plots/PCA_by_SMARCA2_expression.png\n")

# 运行DESeq2
dds <- DESeq(dds)

# 提取结果
res <- results(dds)
res <- as.data.frame(res)

# 这里我们需要直接使用行名作为基因符号，因为我们已经处理好了
res <- rownames_to_column(res, var = "symbol")

# 添加Ensembl ID (如果需要的话，但这次我们主要用Symbol)
# (我们将跳过Ensembl ID的转换，因为我们是从Symbol开始的)
res_ordered <- res[order(res$padj), ]

# -- 提取并保存Top相关基因 --
cat("正在提取与SMARCA2表达最相关的基因...\n")

# 筛选掉NA值
res_filtered <- subset(res_ordered, !is.na(padj))

# 正相关
top_pos_genes <- res_filtered %>% 
  filter(log2FoldChange > 0) %>% 
  arrange(padj, desc(log2FoldChange)) %>% 
  head(15) %>%
  mutate(Regulation = "Positive")

# 负相关
top_neg_genes <- res_filtered %>%
  filter(log2FoldChange < 0) %>%
  arrange(padj, log2FoldChange) %>%
  head(15) %>%
  mutate(Regulation = "Negative")

# 合并
top_genes <- rbind(top_pos_genes, top_neg_genes)

# 打印到控制台
cat("\n--- Top 15 Positively Correlated Genes ---\n")
print(top_pos_genes[, c("symbol", "log2FoldChange", "padj")])
cat("\n--- Top 15 Negatively Correlated Genes ---\n")
print(top_neg_genes[, c("symbol", "log2FoldChange", "padj")])

# 保存到文件
write.csv(top_genes, "results/top_correlated_genes.csv", row.names = FALSE)
cat("\nTop相关基因列表已保存到 r_analysis/results/top_correlated_genes.csv\n")

# 保存完整结果
write.csv(res_ordered, file = "results/all_genes_DESeq2_results.csv", row.names = FALSE)

# 筛选并保存显著结果
significant_res <- subset(res_ordered, padj < 0.05)
write.csv(significant_res, file = "results/significant_genes_padj_0.05.csv", row.names = FALSE)

cat(paste("分析完成，发现", nrow(significant_res), "个与SMARCA2表达显著相关的基因。\n"))
cat("完整结果已保存到 r_analysis/results/all_genes_DESeq2_results.csv\n")
cat("显著基因结果已保存到 r_analysis/results/significant_genes_padj_0.05.csv\n")

# ---- 5. 结果可视化 ----

# 准备用于绘图的数据
plot_data <- as.data.frame(res)
plot_data$significant <- ifelse(!is.na(plot_data$padj) & plot_data$padj < 0.05, "Significant", "Not Significant")
plot_data$symbol <- rownames(res)

# 火山图 (出版级质量)
volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  geom_text_repel(
    data = top_genes, # 只标注top基因
    aes(label = gsub("\\|.*$", "", symbol)), # 只显示基因名，不显示ID
    size = 4,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = 20 # 增加可重叠的标签数
  ) +
  scale_color_manual(name = "", values = c("Not Significant" = "grey80", "Significant" = "#c0392b")) +
  theme_bw(base_size = 16) +
  labs(
    title = "与SMARCA2表达相关的基因",
    subtitle = "火山图",
    x = bquote(~Log[2]~ "Fold Change"),
    y = bquote(~-Log[10]~ "(Adjusted P-value)")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey30") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

ggsave("plots/volcano_plot.png", plot = volcano_plot, width = 12, height = 10, dpi = 300)
cat("出版级质量的火山图已保存到 r_analysis/plots/volcano_plot.png\n")

# 提取用于热图的基因
if (nrow(significant_res) > 0) {
    # vsd 已经从PCA部分计算得出，无需重复计算
    top_genes_for_heatmap <- top_genes$symbol
    
    # 确保这些基因在vsd对象的行名中
    top_genes_for_heatmap <- intersect(top_genes_for_heatmap, rownames(assay(vsd)))
    
    if (length(top_genes_for_heatmap) > 1) {
        heatmap_matrix <- assay(vsd)[top_genes_for_heatmap, ]
        
        # 清理行名用于显示
        rownames(heatmap_matrix) <- gsub("\\|.*$", "", rownames(heatmap_matrix))

        # 准备列注释 (样本信息)
        annotation_col <- data.frame(
            SMARCA2_Expression = phenotype_data$SMARCA2_expression
        )
        rownames(annotation_col) <- rownames(phenotype_data)
        
        # 定义热图颜色
        heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
        
        # 为了避免极端值影响整体的可视化效果，我们为颜色设定一个对称的范围
        # Z-score在此范围之外的值将被视为与范围边界值相同
        color_break_limit <- 2
        heatmap_breaks <- seq(-color_break_limit, color_break_limit, length.out = length(heatmap_colors) + 1)
        
        # 绘制热图 (出版级质量)
        pheatmap(heatmap_matrix,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 show_colnames = FALSE,
                 scale = "row", # 按行标准化，以观察相对表达变化
                 annotation_col = annotation_col,
                 color = heatmap_colors,
                 breaks = heatmap_breaks, # <-- 关键修改：手动设置颜色断点
                 main = paste("Top", length(top_genes_for_heatmap), "个与SMARCA2表达最相关的基因"),
                 fontsize = 12,
                 fontsize_row = 10,
                 filename = "plots/top_genes_heatmap.png",
                 width = 12,
                 height = 10)
                 
        cat(paste("Top", length(top_genes_for_heatmap), "个基因的热图已保存到 r_analysis/plots/top_genes_heatmap.png\n"))
    } else {
        cat("没有足够的Top基因来绘制热图。\n")
    }
} else {
    cat("没有足够的显著基因来绘制热图。\n")
}

# ---- 6. 功能富集分析 (GO & KEGG) ----
if (nrow(significant_res) > 0) {
    cat("开始对显著基因进行功能富集分析...\n")
    
    # 我们需要从原始基因ID中提取出Gene Symbol或EntrezID
    # 这里我们提取EntrezID，因为它更稳定
    entrez_ids_raw <- gsub(".*\\|", "", significant_res$symbol)
    entrez_ids <- entrez_ids_raw[!is.na(as.numeric(entrez_ids_raw))] # 确保是数字ID
    
    # GO 富集分析
    go_results <- enrichGO(gene          = entrez_ids,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENTREZID',
                           ont           = "ALL", # "BP", "MF", "CC"
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.2)

    if (!is.null(go_results) && nrow(go_results) > 0) {
        # 保存GO结果
        write.csv(as.data.frame(go_results), "results/GO_enrichment_results.csv", row.names = FALSE)
        cat("GO富集分析结果已保存到 r_analysis/results/GO_enrichment_results.csv\n")
        
        # 可视化GO结果
        go_plot <- dotplot(go_results, showCategory = 30) + labs(title = "GO 富集分析")
        ggsave("plots/GO_enrichment_dotplot.png", plot = go_plot, width = 12, height = 10)
        cat("GO富集分析点状图已保存到 r_analysis/plots/GO_enrichment_dotplot.png\n")
    } else {
        cat("GO富集分析未找到显著结果。\n")
    }
    
    # KEGG 通路富集分析
    kegg_results <- enrichKEGG(gene         = entrez_ids,
                               organism     = 'hsa', # hsa for human
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05)

    if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
        # 保存KEGG结果
        write.csv(as.data.frame(kegg_results), "results/KEGG_enrichment_results.csv", row.names = FALSE)
        cat("KEGG富集分析结果已保存到 r_analysis/results/KEGG_enrichment_results.csv\n")

        # 可视化KEGG结果
        kegg_plot <- dotplot(kegg_results, showCategory = 30) + labs(title = "KEGG 通路富集分析")
        ggsave("plots/KEGG_enrichment_dotplot.png", plot = kegg_plot, width = 12, height = 10)
        cat("KEGG富集分析点状图已保存到 r_analysis/plots/KEGG_enrichment_dotplot.png\n")
    } else {
        cat("KEGG富集分析未找到显著结果。\n")
    }
} else {
    cat("没有显著基因，跳过功能富集分析。\n")
}

cat("R脚本执行完毕。\n") 