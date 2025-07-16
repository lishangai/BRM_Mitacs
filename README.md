# BRM基因RNA-seq分析项目

## 📋 项目概述

本项目对BRM基因（SMARCA2）进行全面的RNA-seq分析，包括差异表达分析、相关性分析、单细胞RNA-seq分析和生存分析。基于TCGA卵巢癌数据集和单细胞RNA-seq数据，深入探究BRM基因在癌症中的作用机制。

## 🎯 主要功能

- **差异表达分析**: 基于BRM基因表达水平进行分组分析
- **相关性分析**: 识别与BRM基因表达相关的基因
- **单细胞RNA-seq分析**: 细胞质控、聚类、差异表达分析
- **生存分析**: 批量生存分析，识别预后相关基因
- **通路富集分析**: GO和KEGG通路富集分析
- **可视化**: 火山图、热图、生存曲线等多种图表

## 📁 项目结构

```
BRM_Mitacs/
├── data/                           # 数据目录
│   ├── raw/                        # 原始数据
│   │   ├── illuminahiseq_rnaseqv2-RSEM_genes_normalized.csv
│   │   ├── OV.clin.merged.picked.txt
│   │   └── singlecell/
│   │       ├── annotation_HGSC.tsv
│   │       └── UMIcounts_HGSC.tsv
│   └── processed/                  # 处理后数据
│       ├── HGSC_data.h5ad
│       ├── working.csv
│       └── OV.clin.merged.picked.txt
├── src/                           # 分析脚本
│   ├── brm_analysis_csv.py        # 主要RNA-seq分析脚本
│   ├── brm_analysis_pipeline.py   # 单细胞RNA-seq分析脚本
│   ├── batch_survival_analysis.py # 批量生存分析脚本
│   ├── go_enrichment.py          # 通路富集分析脚本
│   ├── visualize_results.py       # 结果可视化脚本
│   └── convert_to_h5ad.py        # 数据格式转换脚本
├── notebooks/                     # Jupyter笔记本
│   ├── BRM基因RNA-seq分析完整流程.ipynb
│   └── BRM_scRNAseq_Analysis.ipynb
├── results/                       # 分析结果
│   ├── brm_analysis_csv/         # RNA-seq分析结果
│   ├── brm_analysis_pipeline/    # 单细胞分析结果
│   ├── batch_survival_analysis/  # 生存分析结果
│   ├── go_enrichment/           # 通路富集结果
│   └── visualize_results/        # 可视化结果
├── docs/                         # 文档
│   ├── BRM_Gene_RNAseq_Analysis_Report.md
│   ├── BRM基因RNA-seq分析综合报告.md
│   └── 质量控制详细对比.md
├── r_analysis/                   # R分析脚本(目前使用不到)
├── environment.yml               # Conda环境配置（推荐）
├── requirements.txt              # pip依赖配置
└── README.md                    # 项目说明
```

## 📥 数据获取

### 方法一：Google Drive下载（推荐）

1. **下载数据文件**：
   - 访问Google Drive链接：[数据文件下载]
   https://drive.google.com/file/d/1HSdLA1D9exaOZ_iLUzfXoyufyCEq8Cme/view?usp=drive_link
   - 下载 `data.zip` 文件

2. **解压数据**：
   ```bash
   # 解压到项目根目录
   unzip data.zip -d .
   ```

3. **验证数据**：
   ```bash
   # 检查数据文件是否存在
   ls data/raw/
   ls data/processed/
   ```

## 🚀 快速开始

### 1. 环境配置

#### 方法一：使用Conda环境（推荐）
```bash
# 创建conda环境
conda env create -f environment.yml

# 激活环境
conda activate brm_analysis
```

#### 方法二：使用pip安装
```bash
# 创建虚拟环境
python -m venv brm_env
source brm_env/bin/activate  # Linux/Mac
# 或
brm_env\Scripts\activate     # Windows

# 安装依赖
pip install -r requirements.txt
```

#### 方法三：手动安装核心依赖
```bash
# 安装基础包
pip install numpy pandas scipy matplotlib seaborn scikit-learn

# 安装分析专用包
pip install scanpy gseapy lifelines statsmodels umap-learn harmonypy

# 安装Jupyter环境
pip install jupyter ipykernel
```

### 2. 验证环境安装
```bash
# 验证Python版本
python --version

# 验证核心包安装
python -c "import numpy, pandas, scanpy, gseapy, lifelines; print('环境配置成功！')"

# 启动Jupyter Lab
jupyter lab
```

### 3. 运行分析

#### RNA-seq差异表达分析
```bash
# 运行主要分析
python src/brm_analysis_csv.py

# 查看结果
ls results/brm_analysis_csv/
```

#### 通路富集分析
```bash
# 运行通路富集分析
python src/go_enrichment.py

# 查看结果
ls results/go_enrichment/
```

#### 生存分析
```bash
# 运行批量生存分析
python src/batch_survival_analysis.py

# 查看结果
ls results/batch_survival_analysis/
```

#### 单细胞RNA-seq分析
```bash
# 运行单细胞分析
python src/brm_analysis_pipeline.py

# 查看结果
ls results/brm_analysis_pipeline/
```


## 📊 分析流程

### RNA-seq分析流程
1. **数据加载与预处理**
   - 加载RSEM标准化数据
   - 过滤低表达基因
   - 数据质量检查

2. **BRM基因分析**
   - 识别BRM/SMARCA2基因
   - 基于表达水平分组（高/低表达）

3. **差异表达分析**
   - Mann-Whitney U检验
   - FDR校正
   - 筛选显著差异基因

4. **相关性分析**
   - Spearman等级相关
   - 识别与BRM相关的基因

### 单细胞RNA-seq分析流程
1. **数据加载与质量控制**
   - 加载H5AD格式数据
   - 细胞和基因过滤
   - 线粒体基因检测

2. **数据预处理**
   - 标准化
   - 高变基因选择
   - PCA降维

3. **批次校正与聚类**
   - Harmony批次校正
   - UMAP可视化
   - K-means聚类

4. **差异表达分析**
   - 细胞类型注释
   - 差异基因识别
   - 通路富集分析

### 生存分析流程
1. **数据准备**
   - 基因表达数据
   - 临床生存数据
   - 样本匹配

2. **批量生存分析**
   - Log-rank检验
   - Cox回归分析
   - 年龄调整分析

3. **结果可视化**
   - Kaplan-Meier曲线
   - 统计报告

## 🔧 环境要求

### Python版本
- Python 3.8-3.11 (推荐3.11)

### 核心依赖包
```
# 基础数据分析
numpy>=1.26.0
pandas>=2.2.0
scipy>=1.15.0
matplotlib>=3.10.0
seaborn>=0.13.0
scikit-learn>=1.7.0

# 单细胞分析
anndata==0.11.4
scanpy==1.11.3

# 通路富集分析
gseapy==1.1.9

# 生存分析
lifelines==0.30.0
statsmodels==0.14.5

# 降维和可视化
umap-learn==0.5.9.post2
harmonypy==0.0.10

# 深度学习
torch==2.7.1

# 数据处理
h5py==3.14.0
tqdm==4.67.1

# Jupyter环境
jupyter>=1.0.0
ipykernel>=6.0.0
```

## 📋 输出结果

### RNA-seq分析结果
- `results/brm_analysis_csv/`
  - `differential_expression_results.csv`: 差异表达结果
  - `correlation_analysis_results.csv`: 相关性分析结果
  - `significant_DEGs.csv`: 显著差异基因
  - `significant_correlations.csv`: 显著相关基因

### 单细胞分析结果
- `results/brm_analysis_pipeline/`
  - `umap_overview.png`: UMAP聚类图
  - `cancer_cells_deg_results.csv`: 癌细胞差异表达结果
  - `smarca2_expression_umap.png`: SMARCA2表达分布图

### 生存分析结果
- `results/batch_survival_analysis/`
  - `all_genes_survival_results.csv`: 所有基因生存分析结果
  - `top_genes_km_curves.png`: Top基因生存曲线
  - `cox_regression_forest_plot.png`: Cox回归森林图

### 通路富集结果
- `results/go_enrichment/`
  - `GO_Biological_Process_results.csv`: GO生物过程结果
  - `KEGG_pathway_results.csv`: KEGG通路结果
  - `enrichment_dotplots.png`: 富集点图

## 📖 详细文档

- **分析报告**: `docs/BRM_Gene_RNAseq_Analysis_Report.md`
- **中文报告**: `docs/BRM基因RNA-seq分析综合报告.md`
- **质量控制**: `docs/质量控制详细对比.md`
- **使用指南**: `docs/guide.md`

## 🔧 故障排除

### 常见环境问题

#### 1. Conda环境创建失败
```bash
# 清理conda缓存
conda clean --all

# 重新创建环境
conda env create -f environment.yml --force
```

#### 2. pip安装包冲突
```bash
# 升级pip
pip install --upgrade pip

# 使用--no-deps安装特定版本
pip install scanpy==1.11.3 --no-deps
```

#### 3. 内存不足
- 减少并行处理数量
- 使用较小的数据子集进行测试
- 增加系统虚拟内存

#### 4. 包版本冲突
```bash
# 检查包版本
pip list | grep scanpy

# 重新安装特定版本
pip uninstall scanpy
pip install scanpy==1.11.3
```


## 🙏 致谢

感谢所有为这个项目做出贡献的研究人员和开发者。

---

**注意**: 请确保在使用前已正确下载和配置数据文件。如有问题，请参考 `docs/` 目录下的详细文档。