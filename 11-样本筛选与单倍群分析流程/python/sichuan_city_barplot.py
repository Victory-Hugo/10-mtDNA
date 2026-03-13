#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
四川省各市样本数量柱状图生成模块。

功能：
1. 读取样本基本信息 CSV（含 City、Province 列）
2. 筛选四川省（Province == 'Sichuan'）样本
3. 按城市（City）统计样本数量（Sichuan_Count）
4. 绘制柱状图：横轴为 City，纵轴为 Sichuan_Count
5. 输出 PDF 图和可选的统计 CSV

参数化：不硬编码任何路径；兼容 import 调用与命令行调用。
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt


def plot_sichuan_city_barplot(
    input_csv: str,
    output_pdf: str,
    province_col: str = "Province",
    city_col: str = "City",
    province_name: str = "Sichuan",
    output_csv: Optional[str] = None,
    verbose: bool = False,
) -> Tuple[Optional[Path], Path]:
    """绘制四川省各市样本数量柱状图。

    参数：
        input_csv: 样本基本信息 CSV 路径（含 City、Province 列）
        output_pdf: 输出柱状图 PDF 路径
        province_col: 省份列名（默认 'Province'）
        city_col: 城市列名（默认 'City'）
        province_name: 目标省份英文名称（默认 'Sichuan'）
        output_csv: 可选；统计结果 CSV 输出路径
        verbose: 是否打印详细日志

    返回：
        (csv_path_or_None, pdf_path)
    """
    if verbose:
        print("[步骤1] 读取样本数据…")
    df = pd.read_csv(input_csv, encoding="utf-8", low_memory=False)

    if province_col not in df.columns:
        raise KeyError(f"输入数据缺少省份列：{province_col}")
    if city_col not in df.columns:
        raise KeyError(f"输入数据缺少城市列：{city_col}")

    if verbose:
        print(f"  总样本数：{len(df)}")

    if verbose:
        print(f"[步骤2] 筛选 {province_name} 省样本…")
    df_sc = df[df[province_col] == province_name].copy()

    if verbose:
        print(f"  {province_name} 省样本数：{len(df_sc)}")

    if len(df_sc) == 0:
        raise ValueError(
            f"在 {province_col} 列中未找到 '{province_name}' 的记录，"
            f"请确认 province_name 参数是否正确。"
        )

    if verbose:
        print("[步骤3] 按城市统计样本数量…")
    city_counts = (
        df_sc.groupby(city_col, dropna=True)
        .size()
        .reset_index(name="Sichuan_Count")
        .sort_values("Sichuan_Count", ascending=False)
    )

    if verbose:
        print(f"  城市数量：{len(city_counts)}")

    csv_out_path = None
    if output_csv:
        csv_out = Path(output_csv)
        csv_out.parent.mkdir(parents=True, exist_ok=True)
        city_counts.to_csv(csv_out, index=False, encoding="utf-8")
        if verbose:
            print(f"✅ 统计结果保存到：{csv_out}")
        csv_out_path = csv_out

    if verbose:
        print("[步骤4] 绘制柱状图…")

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    fig, ax = plt.subplots(figsize=(max(8, len(city_counts) * 0.6), 6))
    ax.bar(city_counts[city_col], city_counts["Sichuan_Count"], color="#4C72B0")
    ax.set_xlabel("City", fontsize=13)
    ax.set_ylabel("Sichuan_Count", fontsize=13)
    ax.set_title(
        f"Sample counts per city in {province_name} Province", fontsize=14
    )
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    pdf_out = Path(output_pdf)
    pdf_out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(pdf_out)
    plt.close(fig)

    if verbose:
        print(f"✅ 柱状图保存到：{pdf_out}")

    return csv_out_path, pdf_out


def main():
    parser = argparse.ArgumentParser(
        description="绘制四川省各市样本数量柱状图（横轴：City，纵轴：Sichuan_Count）"
    )
    parser.add_argument("--input", "-i", required=True, help="样本基本信息 CSV 路径")
    parser.add_argument(
        "--province-col", "-pc", default="Province", help="省份列名（默认 Province）"
    )
    parser.add_argument(
        "--city-col", "-cc", default="City", help="城市列名（默认 City）"
    )
    parser.add_argument(
        "--province-name",
        "-pn",
        default="Sichuan",
        help="目标省份英文名称（默认 Sichuan）",
    )
    parser.add_argument("--output-pdf", "-op", required=True, help="输出柱状图 PDF 路径")
    parser.add_argument(
        "--output-csv", "-oc", default=None, help="可选；各城市样本数量统计 CSV 输出路径"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="打印详细日志")

    args = parser.parse_args()

    try:
        plot_sichuan_city_barplot(
            input_csv=args.input,
            output_pdf=args.output_pdf,
            province_col=args.province_col,
            city_col=args.city_col,
            province_name=args.province_name,
            output_csv=args.output_csv,
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
