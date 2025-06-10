#!/usr/bin/env python3
"""
测试Excel文件读取
"""

import pandas as pd
import sys

filepath = "illuminahiseq_rnaseqv2-RSEM_genes_normalized (MD5).xlsx"

print("测试文件读取...")
print(f"文件路径: {filepath}")

# 方法1: 默认引擎
try:
    print("\n方法1: pandas默认引擎")
    expr = pd.read_excel(filepath, index_col=0)
    print(f"成功读取! 数据维度: {expr.shape}")
    print(f"前5行前5列:\n{expr.iloc[:5, :5]}")
except Exception as e:
    print(f"失败: {e}")

# 方法2: 指定openpyxl引擎
try:
    print("\n方法2: 指定openpyxl引擎")
    expr = pd.read_excel(filepath, engine='openpyxl', index_col=0)
    print(f"成功读取! 数据维度: {expr.shape}")
    print(f"前5行前5列:\n{expr.iloc[:5, :5]}")
except Exception as e:
    print(f"失败: {e}")

# 方法3: 尝试读取第一个sheet
try:
    print("\n方法3: 读取第一个sheet")
    expr = pd.read_excel(filepath, sheet_name=0, index_col=0)
    print(f"成功读取! 数据维度: {expr.shape}")
    print(f"前5行前5列:\n{expr.iloc[:5, :5]}")
except Exception as e:
    print(f"失败: {e}")

# 方法4: 检查sheet名称
try:
    print("\n方法4: 检查sheet名称")
    xl_file = pd.ExcelFile(filepath)
    print(f"Sheet名称: {xl_file.sheet_names}")
    
    # 读取第一个sheet
    expr = pd.read_excel(xl_file, sheet_name=xl_file.sheet_names[0], index_col=0)
    print(f"成功读取! 数据维度: {expr.shape}")
    print(f"前5行前5列:\n{expr.iloc[:5, :5]}")
except Exception as e:
    print(f"失败: {e}")

print("\n测试完成") 