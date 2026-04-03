#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并模块：将单倍群数据与样本信息合并。

功能：
- 对样本信息进行左连接操作，合并单倍群数据
- 检查缺失的单倍群指派
- 输出合并后的数据

支持两种使用方式：
  1. 命令行调用：python merger.py --haplogroup <csv> --sample <csv> --output <csv>
  2. 模块导入：from merger import merge_haplogroup_with_sample
"""

import argparse
import sys

import pandas as pd


def merge_haplogroup_with_sample(
    haplogroup_csv: str,
    sample_csv: str,
    output_csv: str,
    verbose: bool = False
) -> pd.DataFrame:
    """
    将单倍群数据与样本信息合并。
    
    参数：
        haplogroup_csv: 单倍群CSV文件路径
        sample_csv: 样本信息CSV文件路径
        output_csv: 输出CSV文件路径
        verbose: 是否打印详细日志
    
    返回：
        合并后的DataFrame
    """
    if verbose:
        print("[步骤1] 读取输入文件...")
    
    df_hap = pd.read_csv(haplogroup_csv, encoding='utf-8')
    df_sample = pd.read_csv(sample_csv, encoding='utf-8')
    
    if verbose:
        print(f"  单倍群数据行数：{len(df_hap)}")
        print(f"  样本信息行数：{len(df_sample)}")
    
    if verbose:
        print("[步骤2] 执行左连接（样本信息为左表）...")
    
    df_merged = df_sample.merge(df_hap, how='left', on='ID')
    
    if verbose:
        print(f"  合并后行数：{len(df_merged)}")
    
    if verbose:
        print("[步骤3] 检查缺失的单倍群指派...")
    
    # 检查是否有缺失的单倍群
    haplogroup_cols = [c for c in df_hap.columns if c != 'ID']
    for col in haplogroup_cols:
        if col in df_merged.columns:
            missing_count = df_merged[col].isna().sum()
            if missing_count > 0:
                if verbose:
                    print(f"  ⚠ {col} 中有 {missing_count} 个缺失值")
    
    if verbose:
        print("[步骤4] 去重...")
    
    df_merged = df_merged.drop_duplicates(subset=['ID'])
    
    if verbose:
        print(f"  去重后行数：{len(df_merged)}")
    
    if verbose:
        print("[步骤5] 输出结果...")
    
    df_merged.to_csv(output_csv, index=False, encoding='utf-8')
    
    if verbose:
        print(f"✅ 输出保存到：{output_csv}")
    
    return df_merged


def run(
    haplogroup: str,
    sample: str,
    output: str,
    verbose: bool = False,
) -> int:
    merge_haplogroup_with_sample(
        haplogroup_csv=haplogroup,
        sample_csv=sample,
        output_csv=output,
        verbose=verbose,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    """命令行入口"""
    parser = argparse.ArgumentParser(
        description="合并单倍群数据与样本信息"
    )
    parser.add_argument(
        "--haplogroup",
        required=True,
        help="单倍群CSV文件路径"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="样本信息CSV文件路径"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="输出CSV文件路径"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="打印详细日志"
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        return run(**vars(args))
    except Exception as e:
        print(f"❌ 错误：{e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
