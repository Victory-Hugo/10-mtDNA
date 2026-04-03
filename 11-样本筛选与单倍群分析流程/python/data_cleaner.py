#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
样本筛选和清洗模块：按照原笔记本 (0-样本筛选.ipynb) 的顺序进行操作。

支持两种使用方式：
  1. 命令行调用：python data_cleaner.py --input <xlsx> --mapping <xlsx> --new-source <source> --output <csv>
  2. 模块导入：from data_cleaner import clean_sample_info
"""

import re
import unicodedata
from typing import Optional, List, Dict
import argparse
import sys

import pandas as pd


def standardize_column_from_mapping(
    df: pd.DataFrame,
    col: str,
    mapping_df: pd.DataFrame,
    alias_cols: List[str],
    std_col: str
) -> pd.DataFrame:
    """
    将 df[col] 中的别名替换为映射表中的标准名称。
    
    参数：
        df: 输入DataFrame
        col: 需要标准化的列名
        mapping_df: 映射表DataFrame（已读入内存）
        alias_cols: 别名列名列表
        std_col: 标准名称列名
    
    返回：
        标准化后的DataFrame
    """
    # 构造别名 → 标准名称 映射字典
    alias_to_std = {}
    for alias in alias_cols:
        if alias not in mapping_df.columns:
            raise KeyError(f"映射表中不存在别名列：{alias}")
        alias_to_std.update(zip(mapping_df[alias], mapping_df[std_col]))
    
    # 替换逻辑：缺失值或空字符串保留不动
    def _map_skip_empty(x):
        if pd.isna(x) or (isinstance(x, str) and x.strip() == ''):
            return x
        
        std = alias_to_std.get(x, None)
        if std is None:
            return x
        if pd.isna(std) or (isinstance(std, str) and std.strip() == ''):
            return x
        return std
    
    df_out = df.copy()
    df_out[col] = df_out[col].apply(_map_skip_empty)
    return df_out


def clean_whitespace_and_hidden_chars(df: pd.DataFrame) -> pd.DataFrame:
    """
    清理对象列中的空白和隐形字符。
    
    参数：
        df: 输入DataFrame
    
    返回：
        清理后的DataFrame
    """
    obj_cols = df.select_dtypes(include='object').columns.tolist()
    
    # 隐形字符和控制字符模式
    hidden_ctrl_pattern = (
        r'[\u0000-\u001F\u007F'     # ASCII 控制符
        r'\u200B-\u200D'            # 零宽空格/连字
        r'\u200E\u200F'             # LRM/RLM
        r'\u202A-\u202E'            # 方向控制符
        r'\u2060'                    # word joiner
        r'\uFEFF'                    # BOM
        r'\u00AD'                    # 软连字符
        r']'
    )
    
    def clean_str(s: str) -> str:
        """清理字符串中的空白和隐形字符"""
        # Unicode 规范化
        s = unicodedata.normalize('NFKC', s)
        # NBSP/NNBSP 替换为普通空格
        s = s.replace('\u00A0', ' ').replace('\u202F', ' ')
        # 移除隐形/控制字符
        s = re.sub(hidden_ctrl_pattern, '', s)
        # 移除所有空白字符
        s = re.sub(r'\s+', '', s)
        return s
    
    df_out = df.copy()
    for col in obj_cols:
        s = df_out[col]
        mask = s.notna()
        df_out.loc[mask, col] = s.loc[mask].astype(str).map(clean_str)
    
    return df_out


def clean_sample_info(
    input_excel: str,
    mapping_excel: str,
    new_sample_source: str,
    sample_meta_sheet: str,
    ethnicity_mapping_sheet: str,
    country_mapping_sheet: str,
    province_mapping_sheet: str,
    country_extra_cols: Optional[List[str]] = None,
    ethnicity_extra_cols: Optional[List[str]] = None,
    selected_columns: Optional[List[str]] = None,
    output_csv: Optional[str] = None,
    verbose: bool = False
) -> pd.DataFrame:
    """
    按照原笔记本的顺序进行样本筛选和清洗。
    
    参数：
        input_excel: 输入Excel文件路径
        mapping_excel: 基础信息速查表Excel文件路径
        new_sample_source: 新样本的Source值（用于筛选WHALE样本）
        sample_meta_sheet: 主样本信息 sheet 名
        ethnicity_mapping_sheet: 民族映射 sheet 名
        country_mapping_sheet: 国家映射 sheet 名
        province_mapping_sheet: 省份映射 sheet 名
        country_extra_cols: 从国家映射表中提取的额外列名
        ethnicity_extra_cols: 从民族映射表中提取的额外列名
        selected_columns: 需要选择的列名列表
        output_csv: 输出CSV文件路径
        verbose: 是否打印详细日志
    
    返回：
        清洗后的DataFrame
    """
    
    if verbose:
        print("[步骤1] 读取元数据...")
    
    df_meta = pd.read_excel(input_excel, sheet_name=sample_meta_sheet)
    
    if verbose:
        print(f"  总样本数：{len(df_meta)}")
    
    if verbose:
        print(f"[步骤2] 提取WHALE项目样本（Source == '{new_sample_source}'）...")
    
    df_whale = df_meta.loc[df_meta['Source'] == new_sample_source]
    
    if verbose:
        print(f"  WHALE样本数：{len(df_whale)}")
    
    if verbose:
        print("[步骤3] 提取公共数据（Source不以'LLT'开头）...")
    
    df_public = df_meta.loc[~df_meta['Source'].str.startswith('LLT')]
    df_public = df_public.drop_duplicates(subset=['ID'])
    
    if verbose:
        print(f"  公共样本数（去重后）：{len(df_public)}")
    
    if verbose:
        print("[步骤4] 合并样本并去重...")
    
    df_this = pd.concat([df_whale, df_public], ignore_index=True)
    df_this = df_this.drop_duplicates(subset=['ID'])
    
    if verbose:
        print(f"  合并后样本数：{len(df_this)}")
    
    if verbose:
        print("[步骤5] 读取映射表并标准化处理...")
    
    # 读取映射表
    df_民族 = pd.read_excel(mapping_excel, sheet_name=ethnicity_mapping_sheet)
    df_国家 = pd.read_excel(mapping_excel, sheet_name=country_mapping_sheet)
    df_省份 = pd.read_excel(mapping_excel, sheet_name=province_mapping_sheet)
    
    # 清理映射表的空白
    df_民族 = df_民族.map(lambda x: x.strip() if isinstance(x, str) else x)
    df_国家 = df_国家.map(lambda x: x.strip() if isinstance(x, str) else x)
    df_省份 = df_省份.map(lambda x: x.strip() if isinstance(x, str) else x)
    
    if verbose:
        print("  标准化Population列...")
    
    df_new = standardize_column_from_mapping(
        df_this,
        col='Population',
        mapping_df=df_民族,
        alias_cols=['别名(英文)', '别名(中文)'],
        std_col='标准名称(英文)'
    )
    
    if verbose:
        print("  标准化Country列...")
    
    df_new = standardize_column_from_mapping(
        df_new,
        col='Country',
        mapping_df=df_国家,
        alias_cols=['国家(英文别名)', '国家(中文别名)'],
        std_col='国家(英文标准)'
    )
    
    if verbose:
        print("  标准化Province列...")
    
    df_new = standardize_column_from_mapping(
        df_new,
        col='Province',
        mapping_df=df_省份,
        alias_cols=['省份英文全称(别名)', '省份中文全称(别名)'],
        std_col='省份英文名称(标准)'
    )
    
    if verbose:
        print("[步骤6] 清理空白和隐形字符...")
    
    df_new = clean_whitespace_and_hidden_chars(df_new)
    
    if verbose:
        print("[步骤7] 选择需要的列...")
    
    # 选择列
    if selected_columns is None:
        selected_columns = [
            'ID', 'QC_Other', 'QC_Haplogrep', 'Latitude',
            'Longitude', 'Population', 'City', 'Province', 'Country',
            'Sequence_method', 'Source'
        ]
    
    # 检查列是否存在
    missing_cols = [c for c in selected_columns if c not in df_new.columns]
    if missing_cols:
        if verbose:
            print(f"  警告：以下列不存在，将跳过：{missing_cols}")
        selected_columns = [c for c in selected_columns if c in df_new.columns]
    
    df_new = df_new[selected_columns].copy()
    
    if verbose:
        print("[步骤8] 合并国家和民族信息...")
    
    # 设置默认的额外列
    if country_extra_cols is None:
        country_extra_cols = ['Continent']
    if ethnicity_extra_cols is None:
        ethnicity_extra_cols = ['语系', '语言']
    
    # 确保合并关键列存在
    country_cols_to_select = ['国家(英文标准)'] + country_extra_cols
    country_cols_to_select = [c for c in country_cols_to_select if c in df_国家.columns]
    
    ethnicity_cols_to_select = ['标准名称(英文)'] + ethnicity_extra_cols
    ethnicity_cols_to_select = [c for c in ethnicity_cols_to_select if c in df_民族.columns]
    
    # 合并国家信息
    df_国家_info = df_国家[country_cols_to_select].rename(
        columns={'国家(英文标准)': 'Country'}
    )
    df_new = df_new.merge(df_国家_info, how='left', on='Country').drop_duplicates()
    
    # 合并民族信息
    df_民族_info = df_民族[ethnicity_cols_to_select].rename(
        columns={'标准名称(英文)': 'Population'}
    )
    df_new = df_new.merge(df_民族_info, how='left', on='Population').drop_duplicates()
    df_new = df_new.drop_duplicates(subset=['ID'])
    
    if verbose:
        print("[步骤9] 处理缺失值...")
    
    df_new = df_new.fillna('Unknown')
    df_new.loc[df_new['Province'].isin(['Overseas', 'South']), 'Province'] = 'Unknown'
    
    if verbose:
        print("[步骤10] 重命名语言字段...")
    
    rename_dict = {}
    if '语系' in df_new.columns:
        rename_dict['语系'] = 'Language_family'
    if '语言' in df_new.columns:
        rename_dict['语言'] = 'Language'
    
    if rename_dict:
        df_new = df_new.rename(columns=rename_dict)
    
    if verbose:
        print("[步骤11] 添加分类字段...")
    
    def get_classification(row):
        if row['Country'] == 'China':
            return f"{row['Population']}_{row['Province']}"
        else:
            return f"{row['Population']}_{row['Country']}"
    
    df_new['Classification'] = df_new.apply(get_classification, axis=1)
    
    if verbose:
        print("[步骤12] 输出结果...")
    
    if output_csv:
        df_new.to_csv(output_csv, index=False, encoding='utf-8')
        if verbose:
            print(f"✅ 输出保存到：{output_csv}")
    
    return df_new


def run(
    input: str,
    mapping: str,
    new_source: str,
    sample_meta_sheet: str,
    ethnicity_mapping_sheet: str,
    country_mapping_sheet: str,
    province_mapping_sheet: str,
    country_cols: Optional[List[str]] = None,
    ethnicity_cols: Optional[List[str]] = None,
    output: Optional[str] = None,
    verbose: bool = False,
) -> int:
    clean_sample_info(
        input_excel=input,
        mapping_excel=mapping,
        new_sample_source=new_source,
        sample_meta_sheet=sample_meta_sheet,
        ethnicity_mapping_sheet=ethnicity_mapping_sheet,
        country_mapping_sheet=country_mapping_sheet,
        province_mapping_sheet=province_mapping_sheet,
        country_extra_cols=country_cols,
        ethnicity_extra_cols=ethnicity_cols,
        output_csv=output,
        verbose=verbose,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="样本筛选和清洗：按照原笔记本的顺序进行操作"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入Excel文件路径（包含'现代线粒体DNA汇总' sheet）"
    )
    parser.add_argument(
        "--mapping", "-m",
        required=True,
        help="基础信息速查表Excel文件路径"
    )
    parser.add_argument(
        "--new-source", "-s",
        required=True,
        help="新样本的Source值（用于筛选WHALE样本）"
    )
    parser.add_argument(
        "--sample-meta-sheet",
        required=True,
        help="主样本信息 sheet 名"
    )
    parser.add_argument(
        "--ethnicity-mapping-sheet",
        required=True,
        help="民族映射 sheet 名"
    )
    parser.add_argument(
        "--country-mapping-sheet",
        required=True,
        help="国家映射 sheet 名"
    )
    parser.add_argument(
        "--province-mapping-sheet",
        required=True,
        help="省份映射 sheet 名"
    )
    parser.add_argument(
        "--country-cols", "-cc",
        nargs='+',
        default=['Continent'],
        help="从国家映射表中提取的列名（除了'国家(英文标准)' —— 该列固定）"
    )
    parser.add_argument(
        "--ethnicity-cols", "-ec",
        nargs='+',
        default=['语系', '语言'],
        help="从民族映射表中提取的列名（除了'标准名称(英文)' —— 该列固定）"
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
