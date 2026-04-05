#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
频率分析前置处理模块：准备数据用于频率计算和PCA分析。


支持两种使用方式：
  1. 命令行调用：python frequency_preprocessor.py --input <csv> --output <csv> ...
  2. 模块导入：from frequency_preprocessor import preprocess_for_frequency
"""

import argparse
import sys
from typing import Optional, List

import pandas as pd


def _build_frequency_matrix(df: pd.DataFrame, group_col: str, hap_col: str, min_count: int = 0) -> pd.DataFrame:
    """从数据构建单倍群频率透视表。

    参数：
        df: 输入 DataFrame（已含分组列和单倍群列）
        group_col: 分组列名（如 'Classification'）
        hap_col: 单倍群列名
        min_count: 若 > 0，过滤样本数低于阈值的分组
    返回：
        透视表 DataFrame（行=分组，列=单倍群，值=频率）
    """
    if min_count > 0:
        counts = df.groupby(group_col)[group_col].transform('count')
        df = df[counts >= min_count]
    df_freq = df.groupby([group_col, hap_col]).size().reset_index(name='Count')
    df_freq['Frequency'] = df_freq.groupby(group_col)['Count'].transform(lambda x: x / x.sum())
    df_pivot = df_freq.pivot(
        index=group_col,
        columns=hap_col,
        values='Frequency'
    ).fillna(0).reset_index()
    return df_pivot


def preprocess_for_frequency(
    input_csv: str,
    province_mapping_excel: str,
    selected_columns: Optional[List[str]] = None,
    min_sample_count: int = 20,
    haplogroup_column: str = 'Haplogroup_YuChunLi',
    output_prepared_csv: Optional[str] = None,
    output_prepared_china_csv: Optional[str] = None,
    output_frequency_csv: Optional[str] = None,
    output_frequency_china_csv: Optional[str] = None,
    verbose: bool = False
) -> tuple:
    """
    为频率分析准备数据。
    
    参数：
        input_csv: 输入CSV文件路径（sample_with_haplogroup.csv）
        province_mapping_excel: 省份映射表路径（用于Sinitic分组）
        selected_columns: 需要选择的列名列表
        min_sample_count: Classification样本数最小阈值（默认20）
        haplogroup_column: 单倍群列名（默认'Haplogroup_YuChunLi'）
        output_prepared_csv: 准备好的数据输出路径（全球）
        output_prepared_china_csv: 准备好的数据输出路径（中国）
        output_frequency_csv: 频率透视表输出路径（全球）
        output_frequency_china_csv: 频率透视表输出路径（中国）
        verbose: 是否打印详细日志
    
    返回：
        (df_prepared, df_frequency, df_china_prepared, df_china_frequency)
    """
    
    if verbose:
        print("[步骤1] 读取输入数据...")
    
    df = pd.read_csv(input_csv, encoding='utf-8', low_memory=False)
    
    if verbose:
        print(f"  原始数据行数：{len(df)}")
    
    if verbose:
        print("[步骤2] 选择需要的列...")
    
    # 默认选择的列
    if selected_columns is None:
        selected_columns = [
            'ID', 'Population', 'Province', 'Country',
            haplogroup_column, 'Continent', 'Language',
            'Sequence_method', 'Classification', 'QC_Haplogrep',
            'Language_family'
        ]
    
    # 检查列是否存在
    missing_cols = [c for c in selected_columns if c not in df.columns]
    if missing_cols:
        if verbose:
            print(f"  警告：以下列不存在，将跳过：{missing_cols}")
        selected_columns = [c for c in selected_columns if c in df.columns]
    
    df = df.loc[:, selected_columns]
    
    if verbose:
        print("[步骤3] 替换'Unknown'为NA...")
    
    df = df.replace('Unknown', pd.NA)
    
    if verbose:
        print("[步骤4] 删除Population列为空的样本...")
    
    df = df.dropna(subset=['Population'])
    
    if verbose:
        print(f"  删除后样本数：{len(df)}")
    
    if verbose:
        print("[步骤5] 处理中国和海外数据...")
    
    # 中国数据：必须有Province
    df_china = df.loc[df['Country'] == 'China']
    df_china = df_china.dropna(subset=['Province'])
    
    if verbose:
        print(f"  中国样本数：{len(df_china)}")
    
    # 海外数据：必须有Country
    df_overseas = df.loc[df['Country'] != 'China']
    df_overseas = df_overseas.dropna(subset=['Country'])
    
    if verbose:
        print(f"  海外样本数：{len(df_overseas)}")
    
    # 合并
    df = pd.concat([df_china, df_overseas], ignore_index=True)
    
    if verbose:
        print(f"  合并后样本数：{len(df)}")
    
    if verbose:
        print("[步骤6] 处理Language_family分组（Sinitic north/south）...")
    
    # 读取省份映射表
    df_province_map = pd.read_excel(
        province_mapping_excel,
        sheet_name='省份中华英文'
    ).loc[:, ['省份英文名称(标准)', 'South_North(Qinling-Huaihe)']]
    
    # 合并省份信息
    df = df.merge(
        df_province_map,
        left_on='Province',
        right_on='省份英文名称(标准)',
        how='left'
    ).drop(columns=['省份英文名称(标准)'])
    
    # 创建Group列：Sinitic分为north/south，其他用Language_family
    def create_group(row):
        if pd.isna(row['Language']):
            return pd.NA
        
        if row['Language'] == 'Sinitic':
            if pd.isna(row['South_North(Qinling-Huaihe)']):
                return 'Sinitic'
            elif row['South_North(Qinling-Huaihe)'] == 'South':
                return 'Sinitic_South'
            elif row['South_North(Qinling-Huaihe)'] == 'North':
                return 'Sinitic_North'
            else:
                return 'Sinitic'
        else:
            return row['Language_family']
    
    df['Group'] = df.apply(create_group, axis=1)
    df = df.drop(columns=['South_North(Qinling-Huaihe)'])
    
    if verbose:
        print("[步骤7] 去重...")
    
    df = df.drop_duplicates(subset=['ID'])
    
    if verbose:
        print(f"  去重后样本数：{len(df)}")
    
    if verbose:
        print("[步骤8] 计算Classification样本数...")
    
    df['Count'] = df.groupby('Classification')['Classification'].transform('count')
    
    if verbose:
        print("[步骤9] 排除样本数小于阈值的Classification...")
    
    initial_count = df['Classification'].nunique()
    df = df[df['Count'] >= min_sample_count]
    final_count = df['Classification'].nunique()
    
    if verbose:
        print(f"  Classification数量：{initial_count} -> {final_count}")
        print(f"  样本数：{len(df)}")
    
    if verbose:
        print("[步骤10] 分离全球和中国数据...")
    
    # 分离中国数据
    df_china = df[df['Country'] == 'China'].copy()
    df_global = df.copy()
    
    if verbose:
        print(f"  全球样本数：{len(df_global)}")
        print(f"  中国样本数：{len(df_china)}")
    
    if verbose:
        print("[步骤11] 计算单倍群频率（全球）...")

    df_pivot_global = _build_frequency_matrix(df_global, 'Classification', haplogroup_column)

    if verbose:
        print(f"  全球频率表行数：{len(df_pivot_global)}")
        print(f"  全球单倍群数：{len(df_pivot_global.columns) - 1}")

    if verbose:
        print("[步骤12] 计算单倍群频率（中国）...")

    df_pivot_china = _build_frequency_matrix(df_china, 'Classification', haplogroup_column)
    
    if verbose:
        print(f"  中国频率表行数：{len(df_pivot_china)}")
        print(f"  中国单倍群数：{len(df_pivot_china.columns) - 1}")
    
    if verbose:
        print("[步骤13] 输出结果...")
    
    if output_prepared_csv:
        df_global.to_csv(output_prepared_csv, index=False, encoding='utf-8')
        if verbose:
            print(f"✅ 全球准备数据已保存到：{output_prepared_csv}")
    
    if output_prepared_china_csv:
        df_china.to_csv(output_prepared_china_csv, index=False, encoding='utf-8')
        if verbose:
            print(f"✅ 中国准备数据已保存到：{output_prepared_china_csv}")
    
    if output_frequency_csv:
        df_pivot_global.to_csv(output_frequency_csv, index=False, encoding='utf-8')
        if verbose:
            print(f"✅ 全球频率表已保存到：{output_frequency_csv}")
    
    if output_frequency_china_csv:
        df_pivot_china.to_csv(output_frequency_china_csv, index=False, encoding='utf-8')
        if verbose:
            print(f"✅ 中国频率表已保存到：{output_frequency_china_csv}")
    
    return df_global, df_pivot_global, df_china, df_pivot_china


def run(
    input: str,
    province_mapping: str,
    selected_cols: Optional[List[str]] = None,
    min_count: int = 20,
    haplogroup_col: str = "Haplogroup_YuChunLi",
    output_prepared: Optional[str] = None,
    output_prepared_china: Optional[str] = None,
    output_frequency: Optional[str] = None,
    output_frequency_china: Optional[str] = None,
    verbose: bool = False,
) -> int:
    preprocess_for_frequency(
        input_csv=input,
        province_mapping_excel=province_mapping,
        selected_columns=selected_cols,
        min_sample_count=min_count,
        haplogroup_column=haplogroup_col,
        output_prepared_csv=output_prepared,
        output_prepared_china_csv=output_prepared_china,
        output_frequency_csv=output_frequency,
        output_frequency_china_csv=output_frequency_china,
        verbose=verbose,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    """命令行入口"""
    parser = argparse.ArgumentParser(
        description="频率分析前置处理：准备数据用于频率计算和PCA分析"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入CSV文件路径（sample_with_haplogroup.csv）"
    )
    parser.add_argument(
        "--province-mapping", "-pm",
        required=True,
        help="省份映射表Excel文件路径（用于Sinitic north/south分组）"
    )
    parser.add_argument(
        "--selected-cols", "-sc",
        nargs='+',
        default=[
            'ID', 'Population', 'Province', 'Country',
            'Haplogroup_YuChunLi', 'Continent', 'Language',
            'Sequence_method', 'Classification', 'QC_Haplogrep',
            'Language_family'
        ],
        help="需要选择的列名列表"
    )
    parser.add_argument(
        "--min-count", "-mc",
        type=int,
        default=20,
        help="Classification样本数最小阈值（默认20）"
    )
    parser.add_argument(
        "--haplogroup-col", "-hc",
        default='Haplogroup_YuChunLi',
        help="单倍群列名（默认'Haplogroup_YuChunLi'）"
    )
    parser.add_argument(
        "--output-prepared", "-op",
        required=True,
        help="准备好的数据输出路径（全球）"
    )
    parser.add_argument(
        "--output-prepared-china", "-opc",
        default=None,
        help="准备好的数据输出路径（中国，可选）"
    )
    parser.add_argument(
        "--output-frequency", "-of",
        required=True,
        help="频率透视表输出路径（全球）"
    )
    parser.add_argument(
        "--output-frequency-china", "-ofc",
        default=None,
        help="频率透视表输出路径（中国，可选）"
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
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
