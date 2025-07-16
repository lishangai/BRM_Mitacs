# r_analysis/convert_clinical_to_csv.R

# 这个脚本用于将从cBioPortal或类似来源下载的、制表符分隔的临床数据.txt文件，
# 转换为一个标准的、逗号分隔的.csv文件，以便于阅读和检查。

# 定义输入和输出文件路径
# 假设原始临床数据文件位于项目的根目录
clinical_txt_file <- "OV.clin.merged.txt"
output_csv_file <- "results/clinical_data_readable.csv"

# 创建输出目录（如果尚不存在）
dir.create("results", showWarnings = FALSE)

# 检查输入文件是否存在
if (!file.exists(clinical_txt_file)) {
    stop(paste("错误：找不到输入的临床数据文件 ->", clinical_txt_file, "\n请确认该文件位于项目的根目录下。"))
}

cat(paste("正在读取临床数据文件:", clinical_txt_file, "\n"))

# 读取制表符分隔的 .txt 文件
# TCGA的临床数据文件通常包含以'#'开头的注释行，我们使用 comment.char = '#' 来忽略它们
clinical_data <- read.delim(clinical_txt_file, 
                            header = TRUE, 
                            sep = "\t", 
                            stringsAsFactors = FALSE, 
                            comment.char = "#",
                            check.names = FALSE) # 防止R自动修改列名

# 将数据框写入CSV文件
write.csv(clinical_data, file = output_csv_file, row.names = FALSE, na = "")

cat(paste("成功！已将临床数据转换为CSV格式并保存到:", output_csv_file, "\n")) 