#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: 0-extract_from_excel.py

功能简介:
从Excel文件的默认工作表中提取ID列和Haplogroup列，生成制表符分割的txt文件。
支持多种列名格式（不区分大小写）：
  - ID列: "sample_id", "sample", "id"
  - Haplogroup列: "haplogroup_full", "haplogroup", "hap"

使用方法:
python3 0-extract_from_excel.py \
  --input-excel /path/to/Table_Haplogroup.xlsx \
  --output-txt /path/to/ID_Hap.txt
"""

import argparse
import pandas as pd
import sys


def find_column(columns, patterns):
    """
    从列名列表中查找匹配指定模式之一的列。
    
    Args:
        columns: 列名列表（所有字母转换为小写进行比对）
        patterns: 要匹配的模式列表
    
    Returns:
        匹配的原始列名，如果未找到则返回None
    """
    columns_lower = {col.lower(): col for col in columns}
    for pattern in patterns:
        if pattern.lower() in columns_lower:
            return columns_lower[pattern.lower()]
    return None


def extract_from_excel(input_excel, output_txt):
    """
    从Excel文件中提取ID和Haplogroup列。
    
    Args:
        input_excel: 输入Excel文件路径
        output_txt: 输出txt文件路径
    
    Returns:
        bool: 操作是否成功
    """
    try:
        # 读取Excel文件的第一个工作表
        print(f"正在读取Excel文件: {input_excel}")
        df = pd.read_excel(input_excel, sheet_name=0)
        
        # 获取所有列名
        columns = df.columns.tolist()
        print(f"找到列名: {columns}")
        
        # 查找ID列
        id_patterns = ["sample_id", "sample", "id"]
        id_column = find_column(columns, id_patterns)
        
        if id_column is None:
            print(f"错误: 未找到ID列。支持的列名: {id_patterns}")
            return False
        
        print(f"找到ID列: {id_column}")
        
        # 查找Haplogroup列
        hap_patterns = ["haplogroup_full", "haplogroup", "hap"]
        hap_column = find_column(columns, hap_patterns)
        
        if hap_column is None:
            print(f"错误: 未找到Haplogroup列。支持的列名: {hap_patterns}")
            return False
        
        print(f"找到Haplogroup列: {hap_column}")
        
        # 提取两列数据
        result_df = df[[id_column, hap_column]].copy()
        result_df.columns = ["ID", "Haplogroup"]
        
        # 写入输出文件（制表符分割）
        print(f"正在写入输出文件: {output_txt}")
        result_df.to_csv(output_txt, sep="\t", index=False)
        
        print(f"成功! 提取了 {len(result_df)} 行数据")
        print(f"样本数据预览:")
        print(result_df.head())
        
        return True
        
    except FileNotFoundError:
        print(f"错误: 输入文件不存在: {input_excel}")
        return False
    except Exception as e:
        print(f"错误: {str(e)}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="从Excel文件中提取ID和Haplogroup列"
    )
    parser.add_argument(
        "--input-excel",
        required=True,
        help="输入Excel文件路径"
    )
    parser.add_argument(
        "--output-txt",
        required=True,
        help="输出txt文件路径"
    )
    
    args = parser.parse_args()
    
    success = extract_from_excel(args.input_excel, args.output_txt)
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
