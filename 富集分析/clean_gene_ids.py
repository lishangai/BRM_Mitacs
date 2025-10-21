#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
清洗基因名：将CSV第一列形如 “GENE|12345” 的值替换为 “GENE”。
默认：读取当前脚本同目录下 "significant_DEGs copy.csv"，输出到 "significant_DEGs_no_id.csv"。
用法示例：
    python clean_gene_ids.py                # 使用默认输入/输出
    python clean_gene_ids.py -i "significant_DEGs copy.csv" -o significant_DEGs_no_id.csv
"""

from pathlib import Path
import csv
import argparse


def strip_id(gene: str) -> str:
    if gene is None:
        return gene
    parts = str(gene).split('|', 1)
    return parts[0].strip()


def main():
    parser = argparse.ArgumentParser(description="移除第一列基因名中的'|编号'部分")
    parser.add_argument("-i", "--input", type=str, default=None, help="输入CSV路径（默认：脚本同目录下 'significant_DEGs copy.csv'）")
    parser.add_argument("-o", "--output", type=str, default=None, help="输出CSV路径（默认：脚本同目录下 'significant_DEGs_no_id.csv'）")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    in_path = Path(args.input) if args.input else (script_dir / "significant_DEGs copy.csv")
    out_path = Path(args.output) if args.output else (script_dir / "significant_DEGs_no_id.csv")

    if not in_path.exists():
        print(f"✗ 未找到输入文件: {in_path}")
        return

    print(f"输入: {in_path}")
    print(f"输出: {out_path}")

    total = 0
    changed = 0

    # 读取并处理
    with open(in_path, 'r', encoding='utf-8', errors='replace', newline='') as fin, \
         open(out_path, 'w', encoding='utf-8', newline='') as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)

        for row in reader:
            if not row:
                writer.writerow(row)
                continue
            total += 1
            old = row[0]
            new = strip_id(old)
            if new != old:
                changed += 1
            row[0] = new
            writer.writerow(row)

    print(f"完成：共处理 {total} 行，其中修改 {changed} 行。")
    print("✓ 已生成去ID版本：", out_path)


if __name__ == "__main__":
    main()
