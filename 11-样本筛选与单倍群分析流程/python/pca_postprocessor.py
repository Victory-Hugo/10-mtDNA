#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA后处理模块：合并PCA结果与元数据，生成可视化输入文件。

功能：
1. 读取PCA结果
2. 提取并合并元数据（Classification, Continent/Province, Group等）
3. 创建可视化用列名
4. 输出用于可视化的PCA结果

支持两种使用方式：
  1. 命令行调用：python pca_postprocessor.py --input <csv> --metadata <csv> ...
  2. 模块导入：from pca_postprocessor import postprocess_pca
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List

import pandas as pd


def postprocess_pca(
    pca_result_csv: str,
    prepared_data_csv: str,
    metadata_columns: Optional[List[str]] = None,
    group_column_name: str = 'Group',
    class_big_column: str = 'Group',
    output_csv: Optional[str] = None,
    verbose: bool = False
) -> pd.DataFrame:
    """
    PCA后处理：合并PCA结果与元数据。
    
    参数：
        pca_result_csv: PCA结果CSV文件路径
        prepared_data_csv: 准备好的数据CSV文件路径（用于提取元数据）
        metadata_columns: 需要合并的元数据列名（如['Classification', 'Continent', 'Group']）
        group_column_name: 原始数据中的分组列名
        class_big_column: 用于大分类的列名（会重命名为'Class_big'）
        output_csv: 输出CSV文件路径
        verbose: 是否打印详细日志
    
    返回：
        处理后的DataFrame
    """
    
    if verbose:
        print("[步骤1] 读取PCA结果...")
    
    df_pca = pd.read_csv(pca_result_csv, encoding='utf-8')
    
    if verbose:
        print(f"  PCA结果行数：{len(df_pca)}")
        print(f"  PCA列数：{len(df_pca.columns)}")
    
    if verbose:
        print("[步骤2] 读取准备数据...")
    
    df_prepared = pd.read_csv(prepared_data_csv, encoding='utf-8')
    
    if verbose:
        print(f"  准备数据行数：{len(df_prepared)}")
    
    if verbose:
        print("[步骤3] 提取元数据...")
    
    # 默认元数据列
    if metadata_columns is None:
        metadata_columns = ['Classification', 'Continent', group_column_name]
    
    # 检查列是否存在
    available_cols = ['Classification']  # Classification必须存在
    for col in metadata_columns:
        if col != 'Classification' and col in df_prepared.columns:
            available_cols.append(col)
    
    # 提取去重的元数据
    df_metadata = df_prepared[available_cols].drop_duplicates(subset=['Classification'])
    
    if verbose:
        print(f"  元数据行数：{len(df_metadata)}")
        print(f"  元数据列：{list(df_metadata.columns)}")
    
    if verbose:
        print("[步骤4] 合并PCA结果和元数据...")
    
    df_merged = df_pca.merge(df_metadata, on='Classification', how='left')
    
    if verbose:
        print(f"  合并后行数：{len(df_merged)}")
    
    if verbose:
        print("[步骤5] 创建可视化列...")
    
    # 创建Class_small（使用Classification）
    df_merged['Class_small'] = df_merged['Classification']
    
    # 重命名class_big_column为Class_big
    if class_big_column in df_merged.columns:
        df_merged = df_merged.rename(columns={class_big_column: 'Class_big'})
        if verbose:
            print(f"  重命名 {class_big_column} -> Class_big")
    
    # 重命名Classification为ID
    df_merged = df_merged.rename(columns={'Classification': 'ID'})
    
    if verbose:
        print("  列名重命名：Classification -> ID")
    
    if verbose:
        print("[步骤6] 输出结果...")
    
    if output_csv:
        df_merged.to_csv(output_csv, index=False, encoding='utf-8')
        if verbose:
            print(f"✅ 输出保存到：{output_csv}")
    
    return df_merged


def main():
    """命令行入口"""
    parser = argparse.ArgumentParser(
        description="PCA后处理：合并PCA结果与元数据"
    )
    parser.add_argument(
        "--pca-result", "-pr",
        required=True,
        help="PCA结果CSV文件路径"
    )
    parser.add_argument(
        "--prepared-data", "-pd",
        required=True,
        help="准备好的数据CSV文件路径（用于提取元数据）"
    )
    parser.add_argument(
        "--metadata-cols", "-mc",
        nargs='+',
        default=['Classification', 'Continent', 'Group'],
        help="需要合并的元数据列名"
    )
    parser.add_argument(
        "--group-column", "-gc",
        default='Group',
        help="原始数据中的分组列名"
    )
    parser.add_argument(
        "--class-big-column", "-cbc",
        default='Group',
        help="用于大分类的列名（会重命名为'Class_big'），默认为'Group'"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="输出CSV文件路径"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="打印详细日志"
    )
    
    args = parser.parse_args()
    
    try:
        postprocess_pca(
            pca_result_csv=args.pca_result,
            prepared_data_csv=args.prepared_data,
            metadata_columns=args.metadata_cols,
            group_column_name=args.group_column,
            class_big_column=args.class_big_column,
            output_csv=args.output,
            verbose=args.verbose
        )
        sys.exit(0)
    except Exception as e:
        print(f"❌ 错误：{e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
