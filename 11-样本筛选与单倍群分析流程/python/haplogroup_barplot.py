#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
单倍群频率柱状图生成模块（按 notebook 步骤重构）：

功能（按顺序）：
1. 读取频率矩阵（全球或中国）
2. 合并元数据（从 prepared 数据提取 `merge_key` 与 `group_col`）
3. 删除全为 0 的数值列（保留 `merge_key` 与 `group_col`）
4. 生成柱状图（按 `group_col` 分组，在组之间插入空行以形成间隔）
5. 输出清理后的频率 CSV 与 PDF 图

参数化：不硬编码任何路径；兼容 import 调用与命令行调用。
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt


def _read_color_mapping(color_csv: str) -> dict:
    """读取颜色映射 CSV，兼容无表头的两列格式。

    返回：{haplogroup: color_hex}
    """
    try:
        df = pd.read_csv(color_csv, header=None, names=["Haplogroup", "Color"], encoding="utf-8")
        if df.shape[1] >= 2:
            return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
    except Exception:
        pass

    # 兜底：尝试带表头读取
    df = pd.read_csv(color_csv, encoding="utf-8")
    if df.shape[1] >= 2:
        return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
    raise ValueError(f"颜色文件格式不符合预期：{color_csv}")


def _clean_zero_columns(df: pd.DataFrame, merge_key: str, group_col: str, verbose: bool = False) -> pd.DataFrame:
    """删除所有数值列中全为 0 的列，保留 `merge_key` 与 `group_col`。"""
    numeric_cols = df.select_dtypes(include=["number"]).columns.difference([group_col, merge_key])
    nonzero_cols = [c for c in numeric_cols if df[c].sum() != 0]
    kept_cols = list(df.columns.difference(numeric_cols)) + nonzero_cols
    df_cleaned = df[kept_cols]
    if verbose:
        removed = set(numeric_cols) - set(nonzero_cols)
        print(f"  删除全 0 数值列：{sorted(removed)}")
    return df_cleaned


def _plot_stacked_bar_with_gaps(
    df_cleaned: pd.DataFrame,
    merge_key: str,
    group_col: str,
    color_mapping: dict,
    output_pdf: str,
) -> Path:
    """按分组插入空行并绘制堆积柱状图。"""
    data = df_cleaned.copy()

    # 按组插入空行
    grouped_blocks = []
    numeric_cols = data.select_dtypes(include=["number"]).columns.tolist()
    all_cols = data.columns.tolist()

    for grp_name, grp_df in data.groupby(group_col):
        grouped_blocks.append(grp_df)
        empty_dict = {}
        for col in all_cols:
            if col in numeric_cols:
                empty_dict[col] = 0
            elif col == group_col:
                empty_dict[col] = grp_name
            elif col == merge_key:
                empty_dict[col] = ""
            else:
                empty_dict[col] = ""
        grouped_blocks.append(pd.DataFrame([empty_dict]))

    spaced_data = pd.concat(grouped_blocks, ignore_index=True)

    # 准备用于绘图的数据
    plot_data = spaced_data.set_index(merge_key).select_dtypes(include=["number"])  # 仅数值列
    colors = [color_mapping.get(h, "#CCCCCC") for h in plot_data.columns]

    # 绘图
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    # 字体可能因系统而异，不强制指定不可用字体
    # plt.rcParams['font.family'] = 'Arial'

    fig, ax = plt.subplots(figsize=(40, 15))
    plot_data.plot(kind='bar', stacked=True, color=colors, ax=ax)
    ax.set_title('The haplogroup distribution of populations', fontsize=16)
    ax.set_xlabel('', fontsize=14)
    ax.set_ylabel('Haplogroup frequency', fontsize=14)
    ax.legend(title='Haplotypes', bbox_to_anchor=(1, 1), loc='upper left')
    plt.tight_layout()

    out = Path(output_pdf)
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out)
    plt.close(fig)
    return out


def run(
    input_csv: str,
    prepared_csv: str,
    color_csv: str,
    group_col: str = "Group",
    output_csv: Optional[str] = None,
    output_pdf: Optional[str] = None,
    verbose: bool = False,
) -> Tuple[Optional[Path], Optional[Path]]:
    """执行频率柱状图生成流程。

    返回：输出 (csv_path, pdf_path)
    """
    if verbose:
        print("[步骤1] 读取频率矩阵…")
    df = pd.read_csv(input_csv, encoding="utf-8")

    if verbose:
        print("[步骤2] 合并元数据（分组列）…")
    df_meta = pd.read_csv(prepared_csv, encoding="utf-8")
    merge_key = "Classification"
    if merge_key not in df.columns:
        raise KeyError(f"频率矩阵缺少合并键列：{merge_key}")
    if merge_key not in df_meta.columns:
        raise KeyError(f"准备数据缺少合并键列：{merge_key}")
    if group_col not in df_meta.columns:
        raise KeyError(f"准备数据缺少分组列：{group_col}")

    df = df.merge(
        df_meta[[merge_key, group_col]].drop_duplicates(subset=[merge_key]),
        on=merge_key,
        how='left'
    )

    if verbose:
        print("[步骤3] 去重…")
    df = df.drop_duplicates(subset=[merge_key])

    if verbose:
        print("[步骤4] 删除全为 0 的数值列…")
    df_cleaned = _clean_zero_columns(df, merge_key=merge_key, group_col=group_col, verbose=verbose)

    csv_out_path = None
    if output_csv:
        csv_out = Path(output_csv)
        csv_out.parent.mkdir(parents=True, exist_ok=True)
        df_cleaned.to_csv(csv_out, index=False, encoding="utf-8")
        if verbose:
            print(f"✅ 清理后的频率表保存到：{csv_out}")
        csv_out_path = csv_out

    pdf_out_path = None
    if output_pdf:
        if verbose:
            print("[步骤5] 读取颜色映射并绘图…")
        color_map = _read_color_mapping(color_csv)
        pdf_out = _plot_stacked_bar_with_gaps(
            df_cleaned=df_cleaned,
            merge_key=merge_key,
            group_col=group_col,
            color_mapping=color_map,
            output_pdf=output_pdf,
        )
        if verbose:
            print(f"✅ 柱状图保存到：{pdf_out}")
        pdf_out_path = pdf_out

    return csv_out_path, pdf_out_path


def main():
    parser = argparse.ArgumentParser(description="单倍群频率柱状图生成（全球/中国）")
    parser.add_argument("--input", "-i", required=True, help="输入频率矩阵 CSV（global/china）")
    parser.add_argument("--prepared", "-p", required=True, help="准备数据 CSV（用于提供分组列）")
    parser.add_argument("--color", "-c", required=True, help="颜色映射 CSV（两列：haplogroup,color）")
    # 合并键固定为 Classification，不提供命令行参数
    parser.add_argument("--group-col", "-gc", default="Group", help="分组列名（默认 Group）")
    parser.add_argument("--output-csv", "-oc", required=False, help="清理后的频率 CSV 输出路径（可选）")
    parser.add_argument("--output-pdf", "-op", required=True, help="柱状图 PDF 输出路径")
    parser.add_argument("--verbose", "-v", action="store_true", help="打印详细日志")

    args = parser.parse_args()

    try:
        run(
            input_csv=args.input,
            prepared_csv=args.prepared,
            color_csv=args.color,
            group_col=args.group_col,
            output_csv=args.output_csv,
            output_pdf=args.output_pdf,
            verbose=args.verbose,
        )
        sys.exit(0)
    except Exception as e:
        print(f"❌ 错误：{e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
