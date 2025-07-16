# BRM基因RNA-seq分析项目

## 项目概述
本项目对BRM基因（SMARCA2）进行RNA-seq差异表达分析和相关性分析，基于illumina HiSeq平台的RSEM标准化数据。

## 项目结构

```
BRM_Mitacs/
├── data/                           # 数据目录
│   ├── raw/                        # 原始数据
│   │   └── illuminahiseq_rnaseqv2-RSEM_genes_normalized (MD5).csv
│   └── processed/                  # 处理后的数据
├── scripts/                        # 分析脚本
│   ├── brm_analysis_simple.py      # 简化版分析脚本
│   ├── brm_analysis_csv.py         # CSV版本分析脚本（推荐）
│   └── visualize_results.py        # 结果可视化脚本
├── results/                        # 分析结果
│   ├── tables/                     # 表格结果
│   │   ├── differential_expression_results.csv
│   │   ├── significant_DEGs.csv
│   │   ├── correlation_analysis_results.csv
│   │   └── significant_correlations.csv
│   └── figures/                    # 图表结果
│       ├── volcano_plot.png        # 火山图
│       ├── deg_heatmap.png         # 差异基因热图
│       ├── smarca2_expression_distribution.png  # 基因表达分布
│       ├── correlation_analysis.png # 相关性分析图
│       └── summary_statistics.png  # 统计汇总图
├── config/                         # 配置文件
│   └── requirements.txt            # Python依赖包
├── doc/                           # 文档
│   ├── BRM_Gene_RNAseq_Analysis_Report.md
│   ├── BRM_分析报告.md
│   └── 质量控制详细对比.md
├── BRM基因RNA-seq分析完整流程.ipynb  # Jupyter分析流程
├── guide.md                       # 分析指南
└── README.md                      # 项目说明（本文件）
```

## 快速开始

### 1. 环境配置
```bash
# 激活conda环境
conda activate f:\githubclone\BRM_Mitacs\.conda

# 安装依赖包
pip install -r config/requirements.txt
```

### 2. 运行分析
```bash
# 进入脚本目录
cd scripts

# 运行主要分析（推荐使用CSV版本）
python brm_analysis_csv.py

# 运行可视化（如果环境支持）
python visualize_results.py
```

### 3. 查看结果
- **表格结果**: `results/tables/`
- **图表结果**: `results/figures/`
- **分析报告**: `doc/`

## 分析流程

### 1. 数据预处理
- 加载RSEM标准化数据
- 过滤低表达基因（在≥10%样本中表达>1）
- 识别BRM/SMARCA2基因

### 2. 差异表达分析
- 基于BRM基因表达中位数分组（高表达vs低表达）
- Mann-Whitney U检验
- Benjamini-Hochberg FDR校正
- 筛选显著差异基因（padj < 0.05）

### 3. 相关性分析
- Spearman等级相关
- 与BRM基因表达的相关性
- FDR校正和阈值筛选

### 4. 结果可视化
- 火山图（Volcano Plot）
- 差异基因热图
- 表达分布图
- 相关性散点图

## 主要发现

### 数据概览
- **总基因数**: 18,551（过滤后）
- **总样本数**: 307
- **分析基因**: SMARCA2|6595

### 差异表达结果
- **显著差异基因**: 6,002个
- **上调基因**: 3,192个
- **下调基因**: 2,810个

### Top差异基因
**上调基因**:
- C10orf81|79949: log2FC=2.17
- GPR12|2835: log2FC=1.91
- HP|3240: log2FC=1.88

**下调基因**:
- SLC6A10P|386757: log2FC=-2.03
- GAGE12D|100132399: log2FC=-1.89

### 相关性分析
- **显著相关基因**: 614个
- **正相关**: 361个
- **负相关**: 253个

## 文件说明

### 脚本文件
- `brm_analysis_csv.py`: 主要分析脚本，支持CSV格式，兼容性最好
- `brm_analysis_simple.py`: 简化版本，减少依赖
- `visualize_results.py`: 结果可视化脚本

### 数据文件
- `data/raw/`: 原始RSEM标准化数据
- `results/tables/`: 分析结果表格
- `results/figures/`: 生成的图表

### 配置文件
- `config/requirements.txt`: Python包依赖

## 注意事项

1. **环境兼容性**: 如果遇到NumPy 2.x兼容性问题，建议使用`brm_analysis_csv.py`
2. **文件格式**: 支持CSV和Excel格式，CSV格式兼容性更好
3. **内存需求**: 47MB数据文件，建议8GB以上内存
4. **运行时间**: 完整分析约需要5-10分钟

## 更新日志



## 联系信息
如有问题，请查看`doc/`目录下的详细分析报告。 

https://geneformer.readthedocs.io/en/latest/getstarted.html