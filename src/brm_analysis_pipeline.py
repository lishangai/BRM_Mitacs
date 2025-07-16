import os
import sys
import argparse
import warnings
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

# --- 辅助类与函数 ---

class Tee:
    """一个辅助类，用于将标准输出同时重定向到控制台和文件。"""
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()  # 确保立即写入
    def flush(self):
        for f in self.files:
            f.flush()

def check_dependencies():
    """检查可选的依赖包是否存在"""
    dependencies = {'harmonypy': False, 'gseapy': False}
    try:
        import harmonypy
        dependencies['harmonypy'] = True
        print("(+) 'harmonypy' 已安装，将用于批次校正。")
    except ImportError:
        print("(-) 'harmonypy' 未安装，将跳过批次校正。")
    try:
        import gseapy
        dependencies['gseapy'] = True
        print("(+) 'gseapy' 已安装，将进行通路富集分析。")
    except ImportError:
        print("(-) 'gseapy' 未安装，将跳过通路富集分析。")
    return dependencies

def save_plot(fig, filename, output_dir):
    """保存图像并关闭窗口以释放内存"""
    output_path = output_dir / filename
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  - 图表已保存: {output_path}")

# --- 分析步骤函数 ---

def step1_load_and_qc(data_path):
    """步骤一：从H5AD文件加载数据并进行质量控制"""
    print("\n--- 步骤一：加载数据与质量控制 ---")
    
    # 从H5AD加载数据，这非常快速且内存友好
    print(f"  > 正在从高效的H5AD文件加载数据:\n    - {data_path}")
    try:
        adata = sc.read_h5ad(data_path)
    except FileNotFoundError:
        print(f"(-) 错误：找不到文件 {data_path}。")
        print("  > 请先运行 'src/convert_to_h5ad.py' 脚本来生成此文件。")
        sys.exit(1) # 退出程序

    print(f"  > 原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 质量控制
    print("\n  > 正在进行质量控制...")
    
    # 1. 标记线粒体基因
    print("    1. 标记线粒体基因...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    mt_gene_count = adata.var['mt'].sum()
    print(f"       - 检测到 {mt_gene_count} 个线粒体基因")
    
    # 2. 计算质量控制指标
    print("    2. 计算质量控制指标...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # 输出质量控制统计信息
    print("       - 质量控制指标统计:")
    print(f"         * 每个细胞平均表达基因数: {adata.obs['n_genes_by_counts'].mean():.1f}")
    print(f"         * 每个细胞平均总UMI数: {adata.obs['total_counts'].mean():.1f}")
    print(f"         * 线粒体基因表达比例平均值: {adata.obs['pct_counts_mt'].mean():.2f}%")
    print(f"         * 线粒体基因表达比例中位数: {adata.obs['pct_counts_mt'].median():.2f}%")
    
    # 3. 开始过滤过程
    print("\n    3. 开始细胞和基因过滤...")
    original_cells = adata.n_obs
    original_genes = adata.n_vars
    
    # 过滤低质量细胞 (表达基因数 < 200)
    print("       a) 过滤表达基因数过少的细胞 (< 200个基因)...")
    cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    cells_after = adata.n_obs
    removed_cells = cells_before - cells_after
    print(f"          - 移除了 {removed_cells} 个细胞 ({removed_cells/cells_before*100:.1f}%)")
    print(f"          - 剩余 {cells_after} 个细胞")
    
    # 过滤低表达基因 (在少于3个细胞中表达)
    print("       b) 过滤低表达基因 (在少于3个细胞中表达)...")
    genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=3)
    genes_after = adata.n_vars
    removed_genes = genes_before - genes_after
    print(f"          - 移除了 {removed_genes} 个基因 ({removed_genes/genes_before*100:.1f}%)")
    print(f"          - 剩余 {genes_after} 个基因")
    
    # 过滤线粒体基因表达比例过高的细胞 (> 20%)
    print("       c) 过滤线粒体基因表达比例过高的细胞 (> 20%)...")
    cells_before = adata.n_obs
    high_mt_cells = (adata.obs['pct_counts_mt'] >= 20).sum()
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    cells_after = adata.n_obs
    print(f"          - 移除了 {high_mt_cells} 个高线粒体表达细胞 ({high_mt_cells/cells_before*100:.1f}%)")
    print(f"          - 剩余 {cells_after} 个细胞")
    
    # 过滤表达基因数过多的细胞 (> 5000个基因，可能是双细胞)
    print("       d) 过滤表达基因数过多的细胞 (> 5000个基因，疑似双细胞)...")
    cells_before = adata.n_obs
    high_gene_cells = (adata.obs['n_genes_by_counts'] >= 5000).sum()
    adata = adata[adata.obs.n_genes_by_counts < 5000, :]
    cells_after = adata.n_obs
    print(f"          - 移除了 {high_gene_cells} 个高基因表达细胞 ({high_gene_cells/cells_before*100:.1f}%)")
    print(f"          - 剩余 {cells_after} 个细胞")
    
    # 总结过滤结果
    total_cells_removed = original_cells - adata.n_obs
    total_genes_removed = original_genes - adata.n_vars
    print(f"\n  > 质量控制总结:")
    print(f"    - 原始数据: {original_cells} 个细胞, {original_genes} 个基因")
    print(f"    - 过滤后数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    print(f"    - 移除细胞: {total_cells_removed} 个 ({total_cells_removed/original_cells*100:.1f}%)")
    print(f"    - 移除基因: {total_genes_removed} 个 ({total_genes_removed/original_genes*100:.1f}%)")
    
    # 输出过滤后的质量指标
    print(f"\n  > 过滤后质量指标:")
    print(f"    - 每个细胞平均表达基因数: {adata.obs['n_genes_by_counts'].mean():.1f} (范围: {adata.obs['n_genes_by_counts'].min()}-{adata.obs['n_genes_by_counts'].max()})")
    print(f"    - 每个细胞平均总UMI数: {adata.obs['total_counts'].mean():.1f} (范围: {adata.obs['total_counts'].min()}-{adata.obs['total_counts'].max()})")
    print(f"    - 线粒体基因表达比例: {adata.obs['pct_counts_mt'].mean():.2f}% (范围: {adata.obs['pct_counts_mt'].min():.2f}%-{adata.obs['pct_counts_mt'].max():.2f}%)")
    
    return adata

def step2_preprocess_and_integrate(adata, use_harmony, output_dir):
    """步骤二：数据预处理、降维和整合"""
    print("\n--- 步骤二：预处理、降维与整合 ---")
    
    # 标准化和高变基因选择
    print("\n  > 1. 数据标准化...")
    print("     - 执行总计数标准化 (target_sum=1e4)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    print(f"       * 标准化后每个细胞总计数: {adata.X.sum(axis=1).mean():.1f}")
    
    print("     - 执行对数变换 (log1p)...")
    sc.pp.log1p(adata)
    print(f"       * 对数变换后表达值范围: {adata.X.min():.2f} - {adata.X.max():.2f}")
    
    print("\n  > 2. 高变基因识别...")
    genes_before_hvg = adata.n_vars
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    hvg_count = adata.var['highly_variable'].sum()
    print(f"     - 参数设置: min_mean=0.0125, max_mean=3, min_disp=0.5")
    print(f"     - 识别出 {hvg_count} 个高变基因 (占总基因的 {hvg_count/genes_before_hvg*100:.1f}%)")
    
    # 检查是否有足够的高变基因
    if hvg_count == 0:
        print("       * 警告: 未识别到高变基因，放宽参数重新尝试...")
        # 放宽参数重新尝试
        sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=5, min_disp=0.3)
        hvg_count = adata.var['highly_variable'].sum()
        print(f"       * 重新识别出 {hvg_count} 个高变基因")
        
        if hvg_count == 0:
            print("       * 错误: 仍然没有高变基因，使用所有基因进行分析")
            adata.var['highly_variable'] = True
            hvg_count = adata.n_vars
    
    elif hvg_count < 500:
        print(f"       * 警告: 高变基因数量较少 ({hvg_count}个)，可能影响下游分析质量")
    
    # 保存原始数据
    print("     - 保存原始标准化数据到 .raw 属性...")
    adata.raw = adata
    
    # 筛选高变基因
    print("     - 筛选数据仅保留高变基因...")
    adata = adata[:, adata.var.highly_variable]
    print(f"     - 筛选后数据维度: {adata.n_obs} 个细胞 × {adata.n_vars} 个基因")
    
    # PCA
    print("\n  > 3. 主成分分析 (PCA)...")
    print("     - 执行数据缩放 (max_value=10)...")
    print(f"       * 缩放前数据维度: {adata.n_obs} 个细胞 × {adata.n_vars} 个基因")
    
    # 检查数据是否为空
    if adata.n_vars == 0:
        print("       * 错误: 没有基因用于PCA分析")
        raise ValueError("没有基因可用于PCA分析，请检查高变基因选择步骤")
    
    if adata.n_obs == 0:
        print("       * 错误: 没有细胞用于PCA分析")
        raise ValueError("没有细胞可用于PCA分析，请检查质量控制步骤")
    
    sc.pp.scale(adata, max_value=10)
    print(f"       * 缩放后表达值范围: {adata.X.min():.2f} - {adata.X.max():.2f}")
    
    print("     - 运行PCA (svd_solver='arpack')...")
    try:
        sc.tl.pca(adata, svd_solver='arpack')
    except Exception as e:
        print(f"       * PCA计算失败: {e}")
        print("       * 尝试使用默认SVD求解器...")
        sc.tl.pca(adata)
    
    explained_variance_ratio = adata.uns['pca']['variance_ratio']
    cumulative_variance = np.cumsum(explained_variance_ratio)
    
    # 安全地计算主成分数量
    pc50_indices = np.where(cumulative_variance >= 0.5)[0]
    pc80_indices = np.where(cumulative_variance >= 0.8)[0]
    
    pc50_index = pc50_indices[0] + 1 if len(pc50_indices) > 0 else len(explained_variance_ratio)
    pc80_index = pc80_indices[0] + 1 if len(pc80_indices) > 0 else len(explained_variance_ratio)
    
    print(f"       * PCA结果: 计算了 {len(explained_variance_ratio)} 个主成分")
    if len(explained_variance_ratio) >= 10:
        print(f"       * 前10个PC解释的方差: {explained_variance_ratio[:10].sum():.1%}")
    else:
        print(f"       * 前{len(explained_variance_ratio)}个PC解释的方差: {explained_variance_ratio.sum():.1%}")
    print(f"       * 前{pc50_index}个PC解释50%方差, 前{pc80_index}个PC解释80%方差")
    print(f"       * 总累积方差: {cumulative_variance[-1]:.1%}")
    
    # 批次校正
    rep_key = 'X_pca'
    if use_harmony:
        print("\n  > 4. Harmony批次校正...")
        print("     - 检查批次变量...")
        if 'sample' in adata.obs.columns:
            sample_counts = adata.obs['sample'].value_counts()
            print(f"       * 发现 {len(sample_counts)} 个样本:")
            for sample, count in sample_counts.items():
                print(f"         - {sample}: {count} 个细胞")
            
            import harmonypy as hm
            print("     - 运行Harmony算法...")
            harmony_out = hm.run_harmony(adata.obsm['X_pca'], adata.obs, vars_use=['sample'])
            adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
            rep_key = 'X_pca_harmony'
            print(f"       * Harmony校正完成，生成了 {adata.obsm['X_pca_harmony'].shape[1]} 维校正后的PCA表示")
        else:
            print("       * 警告: 未找到 'sample' 列，跳过批次校正")
            rep_key = 'X_pca'
    else:
        print("\n  > 4. 跳过批次校正 (Harmony未安装)")

    # 聚类和UMAP
    print(f"\n  > 5. 构建细胞邻近图...")
    print("     - 参数设置: n_neighbors=10, n_pcs=40")
    print(f"     - 使用表示: {rep_key}")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep=rep_key)
    print("       * 邻近图构建完成")
    
    print("\n  > 6. K-Means聚类...")
    from sklearn.cluster import KMeans
    n_clusters = 20
    print(f"     - 参数设置: n_clusters={n_clusters}, random_state=0")
    kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init='auto')
    # 在批次校正后的PCA表示上进行聚类
    cluster_labels = kmeans.fit_predict(adata.obsm[rep_key])
    adata.obs['kmeans'] = cluster_labels.astype(str)
    
    # 统计聚类结果
    cluster_counts = pd.Series(cluster_labels).value_counts().sort_index()
    print(f"     - 聚类结果: 生成了 {len(cluster_counts)} 个聚类")
    print("       * 各聚类细胞数分布:")
    for cluster_id, count in cluster_counts.items():
        print(f"         - 聚类 {cluster_id}: {count} 个细胞")

    print("\n  > 7. UMAP降维可视化...")
    print("     - 运行UMAP算法...")
    sc.tl.umap(adata)
    umap_coords = adata.obsm['X_umap']
    print(f"       * UMAP完成，生成2D坐标")
    print(f"       * X轴范围: [{umap_coords[:, 0].min():.2f}, {umap_coords[:, 0].max():.2f}]")
    print(f"       * Y轴范围: [{umap_coords[:, 1].min():.2f}, {umap_coords[:, 1].max():.2f}]")
    
    # 可视化 - 检查可用的注释列
    print("\n  > 8. 生成UMAP可视化...")
    available_columns = ['kmeans', 'sample', 'status', 'cell_subtype', 'treatment_phase']
    valid_colors = []
    for col in available_columns:
        if col in adata.obs.columns:
            valid_colors.append(col)
            if col != 'kmeans':  # kmeans我们已经知道了
                unique_values = adata.obs[col].nunique()
                print(f"     - {col}: {unique_values} 种不同值")
        else:
            print(f"     - 警告: 列 '{col}' 不存在，跳过")
    
    print(f"     - 将绘制的注释: {', '.join(valid_colors)}")
    sc.pl.umap(adata, color=valid_colors, show=False)
    
    # scanpy在绘制多图时会自动处理figure对象，我们直接保存当前的全局figure
    plt.savefig(output_dir / 'umap_overview.png', dpi=300, bbox_inches='tight')
    plt.close() # 关闭以释放内存
    print(f"  - 图表已保存: {output_dir / 'umap_overview.png'}")
    
    print(f"\n  > 预处理和整合步骤完成!")
    print(f"    - 最终数据维度: {adata.n_obs} 个细胞 × {adata.n_vars} 个高变基因")
    print(f"    - 原始数据维度: {adata.raw.n_vars} 个总基因")
    print(f"    - 可用表示: {list(adata.obsm.keys())}")

    return adata

def step4_differential_analysis(adata, use_gsea, output_dir):
    """步骤四：差异表达与通路富集分析（根据用户提供的新逻辑）"""
    print("\n--- 步骤四：差异表达与通路富集分析 ---")

    # 1. 根据 'cell_subtype' 筛选出卵巢癌细胞
    target_cell_type = 'Ovarian.cancer.cell'
    print(f"\n  > 1. 筛选目标细胞类型...")
    print(f"     - 目标细胞类型: '{target_cell_type}'")
    
    if 'cell_subtype' not in adata.obs.columns:
        print(f"  > (!) 错误: 注释中未找到 'cell_subtype' 列。跳过此步骤。")
        print(f"      可用的观察列: {list(adata.obs.columns)}")
        return None
    
    # 显示所有细胞类型统计
    cell_type_counts = adata.obs['cell_subtype'].value_counts()
    print(f"     - 发现 {len(cell_type_counts)} 种细胞类型:")
    for cell_type, count in cell_type_counts.items():
        print(f"       * {cell_type}: {count} 个细胞")
        
    cancer_cells = adata[adata.obs['cell_subtype'] == target_cell_type].copy()

    if cancer_cells.n_obs == 0:
        print(f"  > (!) 未在 'cell_subtype' 列中找到 '{target_cell_type}'，跳过差异分析。")
        return None
    
    print(f"     - 成功筛选出 {cancer_cells.n_obs} 个卵巢癌细胞")

    # 2. 从癌细胞中，进一步筛选出 'treatment-naive' 的细胞
    target_phase = 'treatment-naive'
    print(f"\n  > 2. 筛选治疗阶段...")
    print(f"     - 目标治疗阶段: '{target_phase}'")
    
    if 'treatment_phase' not in cancer_cells.obs.columns:
        print(f"  > (!) 错误: 注释中未找到 'treatment_phase' 列。跳过此步骤。")
        print(f"      癌细胞中可用的观察列: {list(cancer_cells.obs.columns)}")
        return None
    
    # 显示治疗阶段统计
    treatment_phase_counts = cancer_cells.obs['treatment_phase'].value_counts()
    print(f"     - 在癌细胞中发现 {len(treatment_phase_counts)} 种治疗阶段:")
    for phase, count in treatment_phase_counts.items():
        print(f"       * {phase}: {count} 个细胞")
        
    naive_cancer_cells = cancer_cells[cancer_cells.obs['treatment_phase'] == target_phase].copy()

    if naive_cancer_cells.n_obs == 0:
        print(f"  > (!) 未找到 '{target_phase}' 的癌细胞，跳过差异分析。")
        return None
    
    print(f"     - 成功筛选出 {naive_cancer_cells.n_obs} 个初治癌细胞")

    # 3. 在 'treatment-naive' 的癌细胞中，比较 'Refractory' vs 'Sensitive'
    print(f"\n  > 3. 准备差异表达分析...")
    comparison_groups = {'Refractory', 'Sensitive'}
    
    if 'status' not in naive_cancer_cells.obs.columns:
        print(f"  > (!) 错误: 注释中未找到 'status' 列。跳过此步骤。")
        print(f"      初治癌细胞中可用的观察列: {list(naive_cancer_cells.obs.columns)}")
        return None
    
    # 显示状态统计
    status_counts = naive_cancer_cells.obs['status'].value_counts()
    print(f"     - 在初治癌细胞中发现 {len(status_counts)} 种状态:")
    for status, count in status_counts.items():
        print(f"       * {status}: {count} 个细胞")
        
    available_statuses = set(naive_cancer_cells.obs['status'].unique())
    if not comparison_groups.issubset(available_statuses):
        print("  > (!) 在初治癌细胞中未同时找到 'Refractory' 和 'Sensitive' 组，无法进行比较。")
        print(f"      找到的状态有: {list(available_statuses)}")
        return None

    # 检查病人分布
    if 'patient_id' in naive_cancer_cells.obs.columns:
        print(f"\n     - 病人分布统计:")
        patient_status = naive_cancer_cells.obs.groupby('patient_id')['status'].first()
        refractory_patients = (patient_status == 'Refractory').sum()
        sensitive_patients = (patient_status == 'Sensitive').sum()
        print(f"       * 耐药病人数: {refractory_patients}")
        print(f"       * 敏感病人数: {sensitive_patients}")

    # 执行差异表达分析
    print(f"\n  > 4. 执行差异表达分析...")
    print(f"     - 比较组: 'Refractory' vs 'Sensitive' (对照组)")
    print(f"     - 统计方法: Wilcoxon rank-sum test")
    print(f"     - 分析基因数: {naive_cancer_cells.raw.n_vars}")
    
    # 使用正确的列名 'status' 和组名 'Refractory', 'Sensitive'
    sc.tl.rank_genes_groups(naive_cancer_cells, 'status', groups=['Refractory'], reference='Sensitive', method='wilcoxon')
    
    deg_results = sc.get.rank_genes_groups_df(naive_cancer_cells, group='Refractory')
    
    # 分析差异表达结果
    print(f"\n     - 差异表达分析完成!")
    print(f"       * 总分析基因数: {len(deg_results)}")
    
    # 统计不同显著性水平的基因
    if 'pvals_adj' in deg_results.columns:
        sig_001 = (deg_results['pvals_adj'] < 0.001).sum()
        sig_01 = (deg_results['pvals_adj'] < 0.01).sum()
        sig_05 = (deg_results['pvals_adj'] < 0.05).sum()
        print(f"       * 显著差异基因 (padj < 0.001): {sig_001}")
        print(f"       * 显著差异基因 (padj < 0.01): {sig_01}")
        print(f"       * 显著差异基因 (padj < 0.05): {sig_05}")
    
    # 统计上调和下调基因
    if 'logfoldchanges' in deg_results.columns:
        upregulated = (deg_results['logfoldchanges'] > 0).sum()
        downregulated = (deg_results['logfoldchanges'] < 0).sum()
        print(f"       * 上调基因: {upregulated}")
        print(f"       * 下调基因: {downregulated}")
        
        # 显示top基因
        top_up = deg_results[deg_results['logfoldchanges'] > 0].head(3)
        top_down = deg_results[deg_results['logfoldchanges'] < 0].head(3)
        print(f"       * 前3个上调基因: {', '.join(top_up['names'].tolist())}")
        print(f"       * 前3个下调基因: {', '.join(top_down['names'].tolist())}")
    
    deg_results.to_csv(output_dir / 'cancer_cells_deg_results.csv', index=False)
    print(f"\n     - 差异表达基因结果已保存至: {output_dir / 'cancer_cells_deg_results.csv'}")

    # GSEA分析
    if use_gsea:
        print(f"\n  > 5. 通路富集分析 (GSEA)...")
        import gseapy as gp
        # gene_list 必须是非空且包含数值
        if 'scores' in deg_results.columns and not deg_results.empty:
            gene_list = deg_results.set_index('names')['scores'].sort_values(ascending=False).dropna()
            if not gene_list.empty:
                print(f"     - 基因排序列表长度: {len(gene_list)}")
                print(f"     - Score范围: {gene_list.min():.3f} - {gene_list.max():.3f}")
                print(f"     - 使用基因集: KEGG_2021_Human, GO_Biological_Process_2021")
                
                gsea_results = gp.prerank(rnk=gene_list, gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'], threads=4)
                enriched_terms = len(gsea_results.res2d)
                gsea_results.res2d.to_csv(output_dir / 'gsea_results.csv')
                print(f"     - GSEA分析完成，发现 {enriched_terms} 个富集通路")
                print(f"     - GSEA结果已保存至: {output_dir / 'gsea_results.csv'}")
            else:
                 print("  > (!) 差异基因列表为空，跳过GSEA。")
        else:
            print("  > (!) 差异分析结果中无 'scores' 列，跳过GSEA。")
    else:
        print(f"\n  > 5. 跳过通路富集分析 (gseapy未安装)")
        
    return naive_cancer_cells

def step5_smarca2_analysis(adata, cancer_cells, output_dir):
    """步骤五：聚焦SMARCA2基因分析"""
    print("\n--- 步骤五：SMARCA2基因深度分析 ---")
    target_gene = 'SMARCA2'
    
    print(f"\n  > 1. 检查目标基因...")
    print(f"     - 目标基因: {target_gene}")
    
    if target_gene not in adata.raw.var_names:
        print(f"  > (!) 目标基因 {target_gene} 不在数据中，跳过此步骤。")
        print(f"      数据中包含 {adata.raw.n_vars} 个基因")
        # 查找相似的基因名
        similar_genes = [gene for gene in adata.raw.var_names if 'SMARCA' in gene.upper()]
        if similar_genes:
            print(f"      发现相似基因: {', '.join(similar_genes[:5])}")
        return
    
    print(f"     - 成功找到目标基因 {target_gene}")
    
    # 计算整体表达统计
    smarca2_expression = adata.raw[:, target_gene].X.toarray().flatten()
    print(f"\n  > 2. 全体细胞中{target_gene}表达统计:")
    print(f"     - 表达细胞数: {(smarca2_expression > 0).sum()} / {len(smarca2_expression)} ({(smarca2_expression > 0).mean()*100:.1f}%)")
    print(f"     - 表达均值: {smarca2_expression.mean():.3f}")
    print(f"     - 表达中位数: {np.median(smarca2_expression):.3f}")
    print(f"     - 表达范围: {smarca2_expression.min():.3f} - {smarca2_expression.max():.3f}")

    # 在所有细胞中绘制SMARCA2表达
    print(f"\n  > 3. 生成全体细胞UMAP表达图...")
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color=target_gene, use_raw=True, ax=ax, show=False, title=f'{target_gene} Expression')
    save_plot(fig, 'smarca2_expression_umap.png', output_dir)
    
    # 在我们关心的初治癌细胞中分析
    if cancer_cells is not None and cancer_cells.n_obs > 0:
        print(f"\n  > 4. 初治癌细胞中{target_gene}表达分析:")
        
        # 检查是否有status分组
        if 'status' in cancer_cells.obs.columns:
            status_groups = cancer_cells.obs['status'].unique()
            print(f"     - 发现状态分组: {', '.join(status_groups)}")
            
            # 计算各组表达统计
            cancer_smarca2 = cancer_cells.raw[:, target_gene].X.toarray().flatten()
            for status in status_groups:
                status_mask = cancer_cells.obs['status'] == status
                status_expression = cancer_smarca2[status_mask]
                expressing_cells = (status_expression > 0).sum()
                total_cells = len(status_expression)
                print(f"     - {status}组:")
                print(f"       * 细胞数: {total_cells}")
                print(f"       * 表达细胞: {expressing_cells} ({expressing_cells/total_cells*100:.1f}%)")
                print(f"       * 表达均值: {status_expression.mean():.3f}")
                print(f"       * 表达中位数: {np.median(status_expression):.3f}")
            
            # 绘制小提琴图
            print(f"     - 生成按状态分组的表达小提琴图...")
            fig, ax = plt.subplots(figsize=(6, 5))
            sc.pl.violin(cancer_cells, keys=target_gene, groupby='status', use_raw=True, ax=ax, show=False)
            save_plot(fig, 'smarca2_expression_in_cancer_cells.png', output_dir)
        else:
            print(f"     - 未找到 'status' 列，无法进行分组分析")
    else:
        print(f"\n  > 4. 未提供癌细胞数据，跳过分组分析")

def step5c_patient_level_analysis(cancer_cells, output_dir):
    """(新) 步骤5c：在病人层面比较SMARCA2的表达差异"""
    print("\n--- 步骤5c：病人层面SMARCA2表达分析 ---")
    from scipy.stats import mannwhitneyu

    target_gene = 'SMARCA2'
    
    print(f"\n  > 1. 数据验证...")
    if cancer_cells is None or cancer_cells.n_obs == 0:
        print("  > (!) 缺少癌细胞数据，跳过此分析。")
        return
        
    if target_gene not in cancer_cells.raw.var_names:
        print(f"  > (!) 目标基因 {target_gene} 不在数据中，跳过此步骤。")
        return
    
    print(f"     - 癌细胞数据: {cancer_cells.n_obs} 个细胞")
    print(f"     - 目标基因: {target_gene}")

    # 检查必要的列
    required_columns = ['patient_id', 'status']
    missing_columns = [col for col in required_columns if col not in cancer_cells.obs.columns]
    if missing_columns:
        print(f"  > (!) 缺少必要的列: {', '.join(missing_columns)}")
        return

    # 1. 创建包含病人ID, 状态和SMARCA2表达值的DataFrame
    print(f"\n  > 2. 构建病人层面数据...")
    patient_df = pd.DataFrame({
        'patient_id': cancer_cells.obs['patient_id'],
        'status': cancer_cells.obs['status'],
        'smarca2_expression': cancer_cells.raw[:, target_gene].X.toarray().flatten()
    })
    
    total_cells = len(patient_df)
    unique_patients = patient_df['patient_id'].nunique()
    print(f"     - 总细胞数: {total_cells}")
    print(f"     - 病人数: {unique_patients}")
    
    # 2. 计算每个病人的SMARCA2平均表达
    print(f"\n  > 3. 计算病人层面表达均值...")
    patient_level_data = patient_df.groupby('patient_id').agg(
        mean_smarca2_expression=('smarca2_expression', 'mean'),
        cell_count=('smarca2_expression', 'count'),
        status=('status', 'first') # 每个病人的status是唯一的
    ).reset_index()

    print(f"     - 成功计算了 {len(patient_level_data)} 个病人的表达均值")
    
    # 检查状态分布
    status_patient_counts = patient_level_data['status'].value_counts()
    print(f"     - 病人状态分布:")
    for status, count in status_patient_counts.items():
        print(f"       * {status}: {count} 个病人")

    # 3. 比较Refractory vs Sensitive组的病人
    print(f"\n  > 4. 统计检验...")
    refractory_exp = patient_level_data[patient_level_data['status'] == 'Refractory']['mean_smarca2_expression']
    sensitive_exp = patient_level_data[patient_level_data['status'] == 'Sensitive']['mean_smarca2_expression']

    if len(refractory_exp) > 0 and len(sensitive_exp) > 0:
        print(f"     - 耐药组病人数量: {len(refractory_exp)}")
        print(f"     - 敏感组病人数量: {len(sensitive_exp)}")
        
        # 计算描述性统计
        print(f"     - 耐药组表达统计:")
        print(f"       * 均值: {refractory_exp.mean():.4f}")
        print(f"       * 中位数: {refractory_exp.median():.4f}")
        print(f"       * 标准差: {refractory_exp.std():.4f}")
        
        print(f"     - 敏感组表达统计:")
        print(f"       * 均值: {sensitive_exp.mean():.4f}")
        print(f"       * 中位数: {sensitive_exp.median():.4f}")
        print(f"       * 标准差: {sensitive_exp.std():.4f}")
        
        # 4. 执行Mann-Whitney U检验
        stat, p_value = mannwhitneyu(refractory_exp, sensitive_exp, alternative='two-sided')
        print(f"\n  > 5. 统计检验结果:")
        print(f"     - 检验方法: Mann-Whitney U test (双侧)")
        print(f"     - 统计量: {stat:.2f}")
        print(f"     - P值: {p_value:.6f}")
        
        if p_value < 0.001:
            significance = "*** (p < 0.001)"
        elif p_value < 0.01:
            significance = "** (p < 0.01)"
        elif p_value < 0.05:
            significance = "* (p < 0.05)"
        else:
            significance = "ns (p >= 0.05)"
        print(f"     - 显著性: {significance}")
        
        # 计算效应量
        from scipy.stats import rankdata
        all_values = np.concatenate([refractory_exp, sensitive_exp])
        ranks = rankdata(all_values)
        u1 = stat
        u2 = len(refractory_exp) * len(sensitive_exp) - u1
        effect_size = 1 - (2 * min(u1, u2)) / (len(refractory_exp) * len(sensitive_exp))
        print(f"     - 效应量 (Common Language Effect Size): {effect_size:.4f}")
        
        # 5. 可视化
        print(f"\n  > 6. 生成可视化...")
        fig, ax = plt.subplots(figsize=(6, 6))
        sns.boxplot(data=patient_level_data, x='status', y='mean_smarca2_expression', ax=ax)
        sns.stripplot(data=patient_level_data, x='status', y='mean_smarca2_expression', ax=ax, color='black', size=5)
        ax.set_title(f'SMARCA2 Mean Expression by Patient Status (p={p_value:.4f})')
        ax.set_xlabel('Patient Status')
        ax.set_ylabel(f'Mean {target_gene} Expression')
        save_plot(fig, 'smarca2_patient_level_boxplot.png', output_dir)
    else:
        print("  > (!) 未能同时找到耐药和敏感组的病人，无法进行比较。")
        if len(refractory_exp) == 0:
            print("      缺少耐药组病人")
        if len(sensitive_exp) == 0:
            print("      缺少敏感组病人")
        
    # 保存病人级别的数据
    print(f"\n  > 7. 保存结果...")
    patient_level_data.to_csv(output_dir / 'smarca2_patient_level_analysis.csv', index=False)
    print(f"     - 病人层面分析结果已保存至: {output_dir / 'smarca2_patient_level_analysis.csv'}")
    
    # 保存详细的统计信息
    if len(refractory_exp) > 0 and len(sensitive_exp) > 0:
        summary_stats = {
            'refractory_patients': len(refractory_exp),
            'sensitive_patients': len(sensitive_exp),
            'refractory_mean': refractory_exp.mean(),
            'sensitive_mean': sensitive_exp.mean(),
            'mannwhitney_statistic': stat,
            'p_value': p_value,
            'effect_size': effect_size,
            'significance': significance
        }
        summary_df = pd.DataFrame([summary_stats])
        summary_df.to_csv(output_dir / 'smarca2_statistical_summary.csv', index=False)
        print(f"     - 统计摘要已保存至: {output_dir / 'smarca2_statistical_summary.csv'}")


def step6_prepare_for_geneformer(cancer_cells, deg_results, output_dir):
    """步骤六：为Geneformer模拟准备文件"""
    print("\n--- 步骤六：为Geneformer模拟准备文件 ---")
    
    if cancer_cells is None or deg_results is None:
        print("  > 缺少癌细胞或差异分析结果，跳过此步骤。")
        return
        
    geneformer_dir = output_dir / 'geneformer_input'
    geneformer_dir.mkdir(exist_ok=True)
        
    deg_results.to_csv(geneformer_dir / 'deg_gene_list.csv', index=False)
    
    # 从初治癌细胞中筛选出敏感的作为基线
    sensitive_cancer_cells = cancer_cells[cancer_cells.obs['status'] == 'Sensitive']
    if sensitive_cancer_cells.n_obs > 0:
        baseline_expression = pd.DataFrame(
            sensitive_cancer_cells.raw.X.mean(axis=0).A1,
            index=sensitive_cancer_cells.raw.var_names,
            columns=['mean_expression']
        )
        baseline_expression.to_csv(geneformer_dir / 'baseline_sensitive_cancer_cells.csv')
    else:
        print("  > (!) 未找到敏感的初治癌细胞，无法生成Geneformer的基线表达文件。")
    
    print(f"  > Geneformer输入文件已保存到: {geneformer_dir}")

def main():
    """主执行函数"""
    # --- 1. 设置路径与参数 ---
    PROJECT_ROOT = Path(__file__).resolve().parent.parent
    
    parser = argparse.ArgumentParser(description="单细胞RNA-seq分析流程")
    parser.add_argument('-i', '--input_file', type=str, 
                        default='data/processed/HGSC_data.h5ad',
                        help='处理后的H5AD输入文件路径（相对于项目根目录）')
    parser.add_argument('-o', '--output_dir', type=str, default=None,
                        help='输出目录路径。默认为 results/<脚本名>/')
    args = parser.parse_args()

    # 构建输入路径
    input_path = PROJECT_ROOT / args.input_file
    
    # 构建输出目录
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        script_name = Path(__file__).stem
        output_dir = PROJECT_ROOT / "results" / script_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- 2. 设置日志记录 ---
    log_file_path = output_dir / 'analysis_log.txt'
    original_stdout = sys.stdout
    log_file = open(log_file_path, 'w', encoding='utf-8') # 明确指定UTF-8编码
    sys.stdout = Tee(original_stdout, log_file)

    try:
        print("=" * 60)
        print("=== BRM基因单细胞RNA-seq分析流程 ===")
        print("=" * 60)
        
        # 显示运行环境信息
        print(f"\n> 运行环境信息:")
        print(f"  - Python版本: {sys.version.split()[0]}")
        print(f"  - 工作目录: {os.getcwd()}")
        print(f"  - 项目根目录: {PROJECT_ROOT}")
        print(f"  - 脚本位置: {__file__}")
        
        # 显示配置信息
        print(f"\n> 分析配置:")
        print(f"  - 输入文件: {input_path}")
        print(f"  - 输出目录: {output_dir}")
        print(f"  - 日志文件: {log_file_path}")
        
        # 检查输入文件
        if input_path.exists():
            file_size = input_path.stat().st_size / (1024**2)  # MB
            print(f"  - 输入文件大小: {file_size:.1f} MB")
        else:
            print(f"  - 警告: 输入文件不存在")
        
        # --- 3. 执行分析 ---
        # 忽略一些常见的警告信息和设置Matplotlib
        print(f"\n> 初始化设置...")
        warnings.filterwarnings('ignore')
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['figure.facecolor'] = 'white'
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        sc.settings.verbosity = 1
        print("  - Matplotlib配置完成")
        print("  - Scanpy详细程度设置为1")
        
        print(f"\n> 检查依赖包...")
        dependencies = check_dependencies()
        
        print(f"\n" + "=" * 60)
        print("开始执行分析流程...")
        print("=" * 60)
        
        # 记录开始时间
        import time
        start_time = time.time()
        
        adata = step1_load_and_qc(input_path)
        adata = step2_preprocess_and_integrate(adata, dependencies['harmonypy'], output_dir)
        # 删除了对 step3 的调用
        cancer_cells_adata = step4_differential_analysis(adata, dependencies['gseapy'], output_dir)
        
        if cancer_cells_adata is not None:
            # 确保deg_results是从正确的对象和组中获取的
            deg_results = sc.get.rank_genes_groups_df(cancer_cells_adata, group='Refractory')
            step5_smarca2_analysis(adata, cancer_cells_adata, output_dir)
            step5c_patient_level_analysis(cancer_cells_adata, output_dir) # 调用新的、正确的函数
            step6_prepare_for_geneformer(cancer_cells_adata, deg_results, output_dir)
        
        # 保存最终结果
        print(f"\n--- 保存最终处理数据 ---")
        final_adata_path = output_dir / 'processed_adata.h5ad'
        print(f"  > 保存最终AnnData对象到: {final_adata_path}")
        adata.write(final_adata_path)
        
        # 计算运行时间
        end_time = time.time()
        runtime = end_time - start_time
        print(f"\n" + "=" * 60)
        print(f"=== 分析流程完成 ===")
        print(f"  - 总运行时间: {runtime/60:.1f} 分钟")
        print(f"  - 最终AnnData保存至: {final_adata_path}")
        print(f"  - 所有结果保存至: {output_dir}")
        print("=" * 60)

    except Exception as e:
        print(f"\n" + "!" * 60)
        print(f"!!! 分析流程发生错误 !!!")
        print(f"错误信息: {e}")
        print("!" * 60)
        import traceback
        traceback.print_exc()
    
    finally:
        # --- 4. 清理 ---
        print(f"\n> 清理和结束:")
        print(f"  - 恢复标准输出")
        print(f"  - 关闭日志文件")
        print(f"  - 完整日志已保存至: {log_file_path}")
        sys.stdout = original_stdout
        log_file.close()

if __name__ == "__main__":
    main() 