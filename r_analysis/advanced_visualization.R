## r_analysis/advanced_visualization.R
# ==============================================================================
# 
# 脚本目标:
#   利用 `brm_analysis.R` 生成的核心结果和原始临床数据，
#   创建一套更深入、更直观、面向发表的分析与可视化图表，
#   以全面回答研究目标。
#
# 依赖数据:
#   - `results/all_genes_DESeq2_results.csv`: DESeq2的全部分析结果。
#   - `OV.clin.merged.txt`: 原始临床数据文件 (位于项目根目录)。
#   - `OV.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt`: 
#     原始基因表达数据 (位于 r_analysis 目录)，用于提取SMARCA2原始表达量。
#
# =================================D=============================================

## ---- 1. 环境设置与包加载 ----

# 设置CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# 检查 BiocManager 是否已安装
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# 定义需要的包
required_packages <- c(
  "ggplot2",
  "dplyr",
  "tibble",
  "survival",   # 用于生存分析的核心包
  "survminer",  # 用于绘制美观的生存曲线
  "viridis"     # 用于提供美观的色谱
)

# 安装/加载包
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

# 创建输出目录
dir.create("plots", showWarnings = FALSE)

cat("高级可视化脚本开始运行...\n")


## ---- 2. 数据加载与预处理 ----

# 定义文件路径
deseq_results_file <- "results/all_genes_DESeq2_results.csv"
clinical_file <- "../OV.clin.merged.txt" # 临床数据在项目根目录
counts_file <- "OV.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"

# ---- 2.1 加载 DESeq2 结果 ----
if (!file.exists(deseq_results_file)) {
    stop("错误: 找不到DESeq2结果文件。请先成功运行 brm_analysis.R")
}
deseq_res <- read.csv(deseq_results_file)

# ---- 2.2 加载并清洗临床数据 ----
if (!file.exists(clinical_file)) {
    stop("错误: 找不到临床数据文件 OV.clin.merged.txt。请将其放置在项目根目录。")
}
clinical_data_raw <- read.delim(clinical_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# 清洗：移除信息量低的列
missing_percentage <- sapply(clinical_data_raw, function(x) sum(is.na(x) | x == '' | x == '[Not Available]' | x == '[Not Applicable]') / nrow(clinical_data_raw))
cols_to_remove_missing <- names(which(missing_percentage > 0.8))
cols_to_remove_constant <- names(which(sapply(clinical_data_raw, function(x) length(unique(x))) == 1))
cols_to_remove <- union(cols_to_remove_missing, cols_to_remove_constant)
clinical_data <- clinical_data_raw[, !names(clinical_data_raw) %in% cols_to_remove]

cat("临床数据加载并清洗完毕。\n")

# ---- 2.3 加载原始Counts数据以获取SMARCA2表达量 ----
# (这是一个简化的加载过程，只为提取SMARCA2)
header <- read.delim(counts_file, nrows = 1, header = FALSE, stringsAsFactors = FALSE)
sample_ids_raw <- as.character(header[1, seq(2, ncol(header), by = 3)])
full_counts_data <- read.delim(counts_file, skip = 1, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
raw_count_columns <- seq(1, ncol(full_counts_data), by = 3)
counts_matrix <- full_counts_data[, raw_count_columns]
colnames(counts_matrix) <- sample_ids_raw

# 提取SMARCA2表达
brm_gene_id_raw <- rownames(counts_matrix)[grepl("^SMARCA2\\|", rownames(counts_matrix))]
if(length(brm_gene_id_raw) == 1) {
    smarca2_expression <- as.numeric(counts_matrix[brm_gene_id_raw, ])
    names(smarca2_expression) <- colnames(counts_matrix)
} else {
    stop("找不到或找到多个SMARCA2基因。")
}
cat("SMARCA2原始表达数据已提取。\n")


## ---- 3. 数据整合 ----
cat("开始整合基因表达与临床数据...\n")

# 准备临床数据中的病人ID (通常是大写的 TCGA-XX-YYYY 格式)
clinical_data$patientId_clean <- toupper(clinical_data$patientId)

# 准备表达数据中的病人ID (从长ID中提取 TCGA-XX-YYYY)
expression_patient_ids <- gsub("^(TCGA\\-[0-9A-Z]{2}\\-[0-9A-Z]{4}).*", "\\1", names(smarca2_expression))
expression_patient_ids <- toupper(expression_patient_ids)

# 创建一个用于合并的数据框
merged_df <- data.frame(
  full_sample_id = names(smarca2_expression),
  patientId_clean = expression_patient_ids,
  smarca2_expression = smarca2_expression
)

# 去重：一个病人可能对应多个测序样本（如原发/复发），这里保留表达量最高的那个
merged_df <- merged_df %>%
  group_by(patientId_clean) %>%
  summarise(
    smarca2_expression = max(smarca2_expression, na.rm = TRUE),
    full_sample_id = first(full_sample_id)
  ) %>%
  ungroup()

# 将临床数据合并进来
final_data <- left_join(merged_df, clinical_data, by = "patientId_clean")

# 清理数据：移除生存时间/状态不明确的样本
final_data <- final_data %>%
  filter(!is.na(OS_MONTHS) & !is.na(OS_STATUS) & OS_MONTHS > 0) %>%
  mutate(OS_STATUS = ifelse(OS_STATUS == 'DECEASED', 1, 0))

cat(paste("数据整合完成，得到", nrow(final_data), "个具有完整生存信息的病人数据。\n"))

## ---- 4. 可视化分析 ----

cat("开始生成高级可视化图表...\n")

# ---- 4.1 SMARCA2 表达分布图 ----
plot_dist <- ggplot(final_data, aes(x = "Ovarian Cancer (TCGA)", y = log2(smarca2_expression + 1))) +
  geom_violin(trim = FALSE, fill = "skyblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_bw(base_size = 16) +
  labs(
    title = "SMARCA2 在TCGA卵巢癌队列中的表达分布",
    x = "",
    y = bquote(~Log[2]~ "(Normalized Counts + 1)")
  )
ggsave("plots/SMARCA2_expression_distribution.png", plot = plot_dist, width = 8, height = 7, dpi = 300)
cat("1/5: SMARCA2表达分布图已保存。\n")

# ---- 4.2 Kaplan-Meier 生存分析 ----
# 根据SMARCA2表达中位数将病人分组
final_data$smarca2_group <- ifelse(final_data$smarca2_expression > median(final_data$smarca2_expression), "High", "Low")
final_data$smarca2_group <- factor(final_data$smarca2_group, levels = c("Low", "High"))

# OS 生存曲线
fit_os <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ smarca2_group, data = final_data)
plot_km_os <- ggsurvplot(
  fit_os,
  data = final_data,
  pval = TRUE,
  risk.table = TRUE,
  title = "SMARCA2表达与总生存期(OS)的关系",
  legend.title = "SMARCA2 Expression",
  legend.labs = c("Low Expression", "High Expression"),
  palette = c("#E7B800", "#2E9FDF"),
  ggtheme = theme_bw(base_size = 14)
)
png("plots/kaplan_meier_survival_OS.png", width = 8, height = 8, units = "in", res = 300)
print(plot_km_os, newpage = FALSE)
dev.off()
cat("2/5: OS生存曲线图已保存。\n")

# PFS 生存曲线
if("PFS_MONTHS" %in% names(final_data) && "PFS_STATUS" %in% names(final_data)) {
  final_data <- final_data %>%
    filter(!is.na(PFS_MONTHS) & !is.na(PFS_STATUS) & PFS_MONTHS >= 0) %>%
    mutate(PFS_STATUS = ifelse(PFS_STATUS == 1, 1, 0))
    
  fit_pfs <- survfit(Surv(PFS_MONTHS, PFS_STATUS) ~ smarca2_group, data = final_data)
  plot_km_pfs <- ggsurvplot(
    fit_pfs,
    data = final_data,
    pval = TRUE,
    risk.table = TRUE,
    title = "SMARCA2表达与无进展生存期(PFS)的关系",
    legend.title = "SMARCA2 Expression",
    legend.labs = c("Low Expression", "High Expression"),
    palette = c("#E7B800", "#2E9FDF"),
    ggtheme = theme_bw(base_size = 14)
  )
  png("plots/kaplan_meier_survival_PFS.png", width = 8, height = 8, units = "in", res = 300)
  print(plot_km_pfs, newpage = FALSE)
  dev.off()
  cat("3/5: PFS生存曲线图已保存。\n")
}

# ---- 4.3 SMARCA2 表达与肿瘤分期的关系 ----
stage_col <- "ajcc_pathologic_tumor_stage"
if(stage_col %in% names(final_data)) {
  stage_data <- final_data %>%
    filter(!is.na(.data[[stage_col]]) & .data[[stage_col]] != "[Not Available]") %>%
    mutate(Stage = gsub("Stage (.*)", "\\1", .data[[stage_col]])) %>%
    filter(Stage %in% c("I", "II", "III", "IV"))
  
  plot_stage <- ggplot(stage_data, aes(x = Stage, y = log2(smarca2_expression + 1))) +
    geom_boxplot(aes(fill = Stage), outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_bw(base_size = 16) +
    labs(
      title = "SMARCA2表达量与肿瘤分期的关系",
      x = "AJCC Pathologic Tumor Stage",
      y = bquote(~Log[2]~ "(Normalized Counts + 1)")
    ) +
    theme(legend.position = "none")
  
  # 加一个统计检验
  if(length(unique(stage_data$Stage)) > 1) {
    stat_test <- compare_means(log2(smarca2_expression + 1) ~ Stage, data = stage_data, method = "anova")
    plot_stage <- plot_stage + stat_compare_means(method = "anova", label.y.npc = "top", label.x.npc = "left")
  }
  
  ggsave("plots/expression_by_tumor_stage.png", plot = plot_stage, width = 9, height = 7, dpi = 300)
  cat("4/5: 肿瘤分期关系图已保存。\n")
}

# ---- 4.4 关键基因关联散点图 ----
# 提取Top1正/负相关基因
top_pos_gene_info <- deseq_res %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% head(1)
top_neg_gene_info <- deseq_res %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% head(1)

# 获取它们的表达量
pos_gene_id <- top_pos_gene_info$symbol
neg_gene_id <- top_neg_gene_info$symbol

if(length(pos_gene_id) == 1 && length(neg_gene_id) == 1) {
  pos_gene_expr <- as.numeric(counts_matrix[pos_gene_id, final_data$full_sample_id])
  neg_gene_expr <- as.numeric(counts_matrix[neg_gene_id, final_data$full_sample_id])
  
  scatter_data <- final_data %>%
    mutate(
      pos_gene_expr = log2(pos_gene_expr + 1),
      neg_gene_expr = log2(neg_gene_expr + 1),
      smarca2_log_expr = log2(smarca2_expression + 1)
    )

  # 正相关图
  plot_scatter_pos <- ggplot(scatter_data, aes(x = smarca2_log_expr, y = pos_gene_expr)) +
    geom_point(alpha = 0.6, color = "darkred") +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("SMARCA2 vs.", gsub("\\|.*$", "", pos_gene_id)),
      subtitle = "Top Positively Correlated Gene",
      x = "SMARCA2 Expression (log2)",
      y = paste(gsub("\\|.*$", "", pos_gene_id), "Expression (log2)")
    )
    
  # 负相关图
  plot_scatter_neg <- ggplot(scatter_data, aes(x = smarca2_log_expr, y = neg_gene_expr)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("SMARCA2 vs.", gsub("\\|.*$", "", neg_gene_id)),
      subtitle = "Top Negatively Correlated Gene",
      x = "SMARCA2 Expression (log2)",
      y = paste(gsub("\\|.*$", "", neg_gene_id), "Expression (log2)")
    )
    
  # 合并图
  combined_scatter <- ggarrange(plot_scatter_pos, plot_scatter_neg, ncol = 2)
  ggsave("plots/correlation_scatter_plots.png", plot = combined_scatter, width = 14, height = 7, dpi = 300)
  cat("5/5: 关键基因关联散点图已保存。\n")
}

cat("高级可视化脚本执行完毕。\n") 