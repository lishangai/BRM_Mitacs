import argparse
import time
from pathlib import Path

import pandas as pd
import scipy as sp
import anndata as ad
from pathlib import Path
import time

def convert_tsv_to_h5ad(counts_file, annotation_file, output_file, chunk_size=1000):
    """
    将巨大的UMI计数TSV文件和注释文件转换为H5AD格式，使用分块读取以优化内存。
    该函数现在会根据注释文件中的细胞来筛选表达矩阵中的细胞。
    """
    try:
        # 1. 读取注释文件，并获取有效的细胞名称列表
        print(f"正在读取注释文件: {annotation_file}")
        annotation_df = pd.read_csv(annotation_file, sep='\t', index_col=0)
        valid_cell_names = annotation_df.index.tolist()
        print(f"成功读取 {len(valid_cell_names)} 个细胞的注释。")

        # 2. 准备分块读取巨大的表达矩阵文件
        print(f"准备分块读取表达矩阵: {counts_file}")
        
        # 首先，读取头信息以获取所有细胞的列名，并确定我们需要哪些列
        with open(counts_file, 'r') as f:
            header = f.readline().strip().split('\t')
        
        # 丢弃第一个元素（通常是 'gene' 或 'Ensembl_ID'）
        all_cell_names_in_matrix = header[1:]
        
        # 找出注释文件中的细胞在表达矩阵头部的索引位置
        # 同时，也需要找出在注释文件和表达矩阵中都存在的细胞
        valid_cell_set = set(valid_cell_names)
        matrix_cell_set = set(all_cell_names_in_matrix)
        
        common_cells = [cell for cell in valid_cell_names if cell in matrix_cell_set]
        
        # 对 common_cells 进行排序，确保与 usecols 加载后的顺序一致
        # 这很重要，因为 pandas 内部可能会按索引顺序加载
        common_cells_indices_map = {name: i for i, name in enumerate(all_cell_names_in_matrix)}
        common_cells.sort(key=lambda name: common_cells_indices_map[name])
        
        if len(common_cells) != len(valid_cell_names):
            print("警告: 注释文件中的一些细胞在表达矩阵中找不到或反之。")
            print(f"注释细胞数: {len(valid_cell_names)}, 表达矩阵细胞数: {len(matrix_cell_set)}, 共同细胞数: {len(common_cells)}")
            # 更新注释DF，只保留共同的细胞
            annotation_df = annotation_df.loc[common_cells]
            print(f"注释文件已更新，仅保留 {len(common_cells)} 个共同细胞。")

        # 获取需要保留的列的索引
        cell_indices_to_keep = [all_cell_names_in_matrix.index(cell) for cell in common_cells]
        
        # 我们需要在使用pd.read_csv时传递的列索引，需要加上基因名那一列（索引0）
        column_indices_for_pd = [0] + [i + 1 for i in cell_indices_to_keep]
        # 创建一个列名列表给新的dataframe
        final_column_names = ['gene'] + common_cells


        all_chunks = []
        genes = []
        
        # 使用pandas的TextFileReader进行分块读取
        # 通过 usecols 参数，我们告诉pandas只读取我们感兴趣的列
        reader = pd.read_csv(
            counts_file, 
            sep='\t', 
            chunksize=chunk_size, 
            index_col=0,
            usecols=column_indices_for_pd
        )

        print("开始逐块读取和处理...")
        start_time = time.time()
        i = 0
        for chunk in reader:
            i += 1
            # 删除了有问题的行: chunk = chunk[common_cells]
            # 因为 usecols 已经保证了我们只读取了正确的列。
            # chunk的列顺序由usecols决定，我们已经对common_cells排序以匹配它
            all_chunks.append(sp.sparse.csc_matrix(chunk.values))
            genes.extend(chunk.index.tolist())
            if i % 10 == 0:
                elapsed = time.time() - start_time
                print(f"  已处理 {i * chunk_size} 行基因... ({elapsed:.2f} 秒)")
        
        print("所有分块读取完毕。")

        # 3. 合并所有分块并创建AnnData对象
        print("正在合并稀疏矩阵...")
        csc_matrix = sp.sparse.vstack(all_chunks, format="csc")
        
        print("正在创建AnnData对象...")
        # 矩阵需要转置，以匹配 '细胞 x 基因' 的标准格式
        # 此时的csc_matrix是 (基因 x 细胞)，注释是 (细胞 x 属性)
        
        # 创建一个正确的 var dataframe
        var_df = pd.DataFrame(index=genes)
        
        # 创建AnnData对象，确保obs的索引是正确的
        adata = ad.AnnData(
            csc_matrix.T,
            obs=annotation_df,
            var=var_df
        )
        print(f"AnnData对象创建成功。形状: {adata.shape[0]} 细胞 × {adata.shape[1]} 基因。")

        # 4. 将AnnData对象写入H5AD文件
        print(f"正在将对象写入H5AD文件: {output_file}")
        output_file.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_file, compression="gzip")
        print("✓ H5AD文件转换成功！")

    except FileNotFoundError as e:
        print(f"错误: {e}")
    except Exception as e:
        print(f"在转换过程中发生未知错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    PROJECT_ROOT = Path(__file__).resolve().parent.parent
    
    parser = argparse.ArgumentParser(description="将大型scRNA-seq TSV文件转换为内存高效的H5AD格式。")
    parser.add_argument(
        '--counts', 
        type=str, 
        default='data/raw/singlecell/UMIcounts_HGSC.tsv',
        help='相对于项目根目录的UMI计数文件路径。'
    )
    parser.add_argument(
        '--annotation', 
        type=str, 
        default='data/raw/singlecell/annotation_HGSC.tsv',
        help='相对于项目根目录的细胞注释文件路径。'
    )
    parser.add_argument(
        '--output', 
        type=str, 
        default='data/processed/HGSC_data.h5ad',
        help='相对于项目根目录的输出H5AD文件路径。'
    )
    args = parser.parse_args()

    # 构建绝对路径
    counts_file = PROJECT_ROOT / args.counts
    annotation_file = PROJECT_ROOT / args.annotation
    output_file = PROJECT_ROOT / args.output

    convert_tsv_to_h5ad(counts_file, annotation_file, output_file)

if __name__ == "__main__":
    main()