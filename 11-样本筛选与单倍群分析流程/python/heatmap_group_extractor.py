#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
频率热图分组文件提取模块

按流程：
1) 从准备数据（frequency_prepared*.csv）读取元数据
2) 可选：提取中国子集（Country == 'China'）
3) 提取指定列（默认：Classification, Group），去重
4) 输出到指定目录（默认由管道传入）

双模式：可 import 的 run() + CLI。
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple, Optional

import pandas as pd


def run(
    input_csv: str,
    group_col: str = "Group",
    output_dir: str = "./output/z-score热图",
    suffix: str = "",
    verbose: bool = False,
) -> Tuple[Path, Optional[Path]]:
    """生成热图用分组映射 CSV。

    固定键列为 'Classification'（不允许外部修改）。

    返回：(global_csv_path, china_csv_path|None)
    """
    if verbose:
        print("[步骤1] 读取准备数据…")
    df = pd.read_csv(input_csv, encoding="utf-8", low_memory=False)

    # 校验列
    key_col = "Classification"
    for col in [key_col, group_col]:
        if col not in df.columns:
            raise KeyError(f"输入中缺少列：{col}")
    # 不再依赖国家列，直接基于输入文件提取映射

    if verbose:
        print("[步骤2] 提取全局映射并去重…")
    df_global = df.loc[:, [key_col, group_col]].drop_duplicates(subset=[key_col])

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 文件名沿用 notebook 约定（避免破坏下游）：
    # 可选后缀（如 _China），用于不同输入文件的命名区分
    sfx = suffix if suffix else ""
    global_name = (
        f"Classification_Group_mapping{sfx}.csv"
        if (key_col == "Classification" and group_col == "Group")
        else f"{key_col}_{group_col}_mapping{sfx}.csv"
    )
    global_path = out_dir / global_name
    df_global.to_csv(global_path, index=False, encoding="utf-8")
    if verbose:
        print(f"✅ 全局映射：{global_path}")

    # 不再内置生成中国子集；若需要，可通过传入不同输入文件及 suffix="_China" 实现
    china_path: Optional[Path] = None

    return global_path, china_path


def main():
    parser = argparse.ArgumentParser(description="提取热图分组映射（可选输出后缀，用于区分中国等子集）")
    parser.add_argument("--input", "-i", required=True, help="准备数据 CSV（如 frequency_prepared.csv）")
    parser.add_argument("--group-col", "-gc", default="Group", help="分组列，默认 Group")
    parser.add_argument("--output-dir", "-o", required=True, help="输出目录（将写入映射 CSV）")
    parser.add_argument("--suffix", "-s", default="", help="可选：输出文件名后缀（例如 _China）")
    parser.add_argument("--verbose", "-v", action="store_true", help="打印详细日志")

    args = parser.parse_args()

    try:
        run(
            input_csv=args.input,
            group_col=args.group_col,
            output_dir=args.output_dir,
            suffix=args.suffix,
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
