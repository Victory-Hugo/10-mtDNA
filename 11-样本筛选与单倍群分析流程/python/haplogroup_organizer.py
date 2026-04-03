#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""兼容当前主流程的单倍群整理入口。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

import pandas as pd

from annotate_haplogroups import build_output_rows, load_tree_index, write_output

log = logging.getLogger(__name__)

OUTPUT_FILENAME = "mtDNA_haplogroup_annotation.tsv"
FINAL_FILENAME = "全部单倍群整理.csv"


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)s] %(message)s",
    )


def extract_raw_haplogroups(input_excel: str) -> pd.DataFrame:
    """从当前流程使用的 Excel 中提取样本 ID 与原始单倍群。"""
    log.info("读取质量控制 sheet")
    df_qc = pd.read_excel(input_excel, sheet_name="质量控制")
    df_qc = df_qc.loc[:, ["SampleID", "Haplogroup"]]

    log.info("读取芯片单倍群 sheet")
    df_chip = pd.read_excel(input_excel, sheet_name="芯片单倍群")
    df_chip = df_chip.loc[:, ["SampleID", "Haplogroup"]]

    df_id_hap = pd.concat([df_qc, df_chip], ignore_index=True)
    df_id_hap = df_id_hap.rename(columns={"SampleID": "ID"})
    df_id_hap = df_id_hap[df_id_hap["ID"].notna() & df_id_hap["Haplogroup"].notna()].copy()
    df_id_hap["Haplogroup"] = df_id_hap["Haplogroup"].map(lambda value: str(value).strip())
    df_id_hap = df_id_hap[df_id_hap["Haplogroup"] != ""]
    return df_id_hap


def build_legacy_output(annotation_path: Path, output_dir: Path) -> Path:
    """将新版完整注释结果压缩成当前老流程继续使用的兼容输出。"""
    df_annotation = pd.read_csv(annotation_path, sep="\t", encoding="utf-8-sig")
    df_legacy = df_annotation[
        [
            "ID",
            "Haplogroup_Standardized",
            "Haplogroup_LLT",
            "Haplogroup_YuChunLi",
        ]
    ].rename(columns={"Haplogroup_Standardized": "Haplogroup"})
    final_path = output_dir / FINAL_FILENAME
    df_legacy.to_csv(final_path, index=False, encoding="utf-8")
    return final_path


def run(
    input: str,
    temp_dir: str,
    data_dir: str,
    tree_index: str,
    log_level: str = "INFO",
) -> int:
    configure_logging(log_level)

    input_path = Path(input)
    temp_path = Path(temp_dir)
    output_path = Path(data_dir)
    tree_index_path = Path(tree_index)

    if not input_path.exists():
        raise FileNotFoundError(f"缺失输入 Excel 文件: {input_path}")
    if not tree_index_path.exists():
        raise FileNotFoundError(f"缺失单倍群树索引文件: {tree_index_path}")

    temp_path.mkdir(parents=True, exist_ok=True)
    output_path.mkdir(parents=True, exist_ok=True)

    log.info("步骤1：提取原始单倍群")
    df_raw = extract_raw_haplogroups(str(input_path))
    id_hap_path = temp_path / "ID_Hap.tsv"
    df_raw.to_csv(id_hap_path, sep="\t", index=False, encoding="utf-8")
    log.info("提取完成：%d 条记录 -> %s", len(df_raw), id_hap_path)

    log.info("步骤2：调用新版单倍群层级划分工具")
    input_rows = list(df_raw[["ID", "Haplogroup"]].itertuples(index=False, name=None))
    tree_data = load_tree_index(tree_index_path)
    output_rows = build_output_rows(input_rows, tree_data)
    annotation_path = temp_path / OUTPUT_FILENAME
    write_output(annotation_path, output_rows)
    log.info("注释完成：%s", annotation_path)

    log.info("步骤3：生成兼容当前流程的输出文件")
    final_path = build_legacy_output(annotation_path, output_path)
    log.info("兼容输出完成：%s", final_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "-i", required=True, help="输入 Excel 文件")
    parser.add_argument("--temp-dir", "-t", required=True, help="临时文件目录")
    parser.add_argument("--data-dir", "-d", required=True, help="兼容输出目录")
    parser.add_argument("--tree-index", required=True, help="单倍群树索引 JSON")
    parser.add_argument("--log-level", default="INFO", help="日志级别")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
