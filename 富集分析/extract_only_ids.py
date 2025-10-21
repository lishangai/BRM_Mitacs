#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
从CSV第一列提取管道符后的ID：
    例：CAPN9|10753 -> 10753
默认输入：脚本同目录下 "significant_DEGs copy.csv"
默认输出：脚本同目录下 "significant_DEGs_only_id.csv"
用法：
    python extract_only_ids.py
    python extract_only_ids.py -i "significant_DEGs copy.csv" -o significant_DEGs_only_id.csv
"""

from pathlib import Path
import csv
import argparse

def get_id(gene: str) -> str:
    if gene is None:
        return gene
    s = str(gene)
    if '|' in s:
        return s.split('|', 1)[1].strip()
    return s.strip()  # 如果没有'|'，原样返回


def main():
    parser = argparse.ArgumentParser(description="提取第一列中'|'后的ID")
    parser.add_argument("-i", "--input", type=str, default=None, help="输入CSV路径（默认：脚本同目录 'significant_DEGs copy.csv'）")
    parser.add_argument("-o", "--output", type=str, default=None, help="输出CSV路径（默认：脚本同目录 'significant_DEGs_only_id.csv'）")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    in_path = Path(args.input) if args.input else (script_dir / "significant_DEGs copy.csv")
    out_path = Path(args.output) if args.output else (script_dir / "significant_DEGs_only_id.csv")

    if not in_path.exists():
        print(f"✗ 未找到输入文件: {in_path}")
        return

    print(f"输入: {in_path}")
    print(f"输出: {out_path}")

    total = 0
    changed = 0

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
            new = get_id(old)
            if new != old:
                changed += 1
            row[0] = new
            writer.writerow(row)

    print(f"完成：共处理 {total} 行，其中修改 {changed} 行。")
    print("✓ 已生成仅ID版本：", out_path)

if __name__ == "__main__":
    main()
