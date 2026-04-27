#!/usr/bin/env python3
"""将 Phylotree Excel 层级表转换为与模板一致的 JSON。

示例：
    python python/1-1-excel_to_phylotree_json.py \
        --input-excel input/线粒体DNA系统发育树.xlsx \
        --template-json template/phylotree_index_withacc.json \
        --output-json output/phylotree_build_18.json
"""

from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Any

import pandas as pd

LOG = logging.getLogger(__name__)


class PhylotreeConvertError(RuntimeError):
    """Phylotree Excel 转换失败。"""


def normalize_text(value: Any) -> str:
    """统一文本空白，去除前后空格。"""
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).replace("\u00a0", " ")).strip()


def first_nonempty_level(row: pd.Series, level_columns: list[str]) -> tuple[int, str] | None:
    """返回一行中第一个非空 Level 单元格的列号与节点名。"""
    for index, column in enumerate(level_columns, start=1):
        value = normalize_text(row.get(column))
        if value:
            return index, value
    return None


def extract_accessions(row: pd.Series, accession_columns: list[str]) -> list[str]:
    """提取 accession 列，保留 Excel 中从左到右的顺序。"""
    accessions: list[str] = []
    for column in accession_columns:
        value = normalize_text(row.get(column))
        if value:
            accessions.append(value)
    return accessions


def build_haplogroups(df: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """从 Phylotree 层级表重建 haplogroups 字典。"""
    level_columns = [f"Level{i}" for i in range(1, 26)]
    accession_columns = [column for column in df.columns if str(column).startswith("Accession ID")]

    missing_level_columns = [column for column in level_columns if column not in df.columns]
    if missing_level_columns:
        raise PhylotreeConvertError(f"输入表缺少 Level 列: {missing_level_columns}")

    level_stack: dict[int, str] = {}
    haplogroups: dict[str, dict[str, Any]] = {}

    for excel_row, (_, row) in enumerate(df.iterrows(), start=2):
        node = first_nonempty_level(row, level_columns)
        if node is None:
            continue

        level_column_number, haplogroup = node
        level = 0 if haplogroup == "mt-MRCA(RSRS)" else level_column_number
        mutation_column = f"Level{level_column_number + 1}"
        mutations = normalize_text(row.get(mutation_column)) if mutation_column in df.columns else ""
        accessions = extract_accessions(row, accession_columns)

        if haplogroup in haplogroups:
            raise PhylotreeConvertError(f"发现重复节点名: {haplogroup} (Excel row {excel_row})")

        if level == 0:
            parent = ""
            lineage = [haplogroup]
        else:
            parent = level_stack.get(level - 1)
            if parent is None:
                raise PhylotreeConvertError(
                    f"无法为节点 {haplogroup} 找到上一级父节点 "
                    f"(Excel row {excel_row}, level {level})"
                )
            lineage = [level_stack[i] for i in range(0, level) if i in level_stack] + [haplogroup]

        record: dict[str, Any] = {}
        if accessions:
            record["accessions"] = accessions
        record["level"] = level
        record["lineage"] = lineage
        record["mutations"] = mutations
        record["parent"] = parent
        haplogroups[haplogroup] = record

        level_stack[level] = haplogroup
        for stacked_level in list(level_stack):
            if stacked_level > level:
                del level_stack[stacked_level]

    if "mt-MRCA(RSRS)" not in haplogroups:
        raise PhylotreeConvertError("未在 Excel 中找到根节点 mt-MRCA(RSRS)")

    return {key: haplogroups[key] for key in sorted(haplogroups)}


def load_template(template_json: Path) -> dict[str, Any]:
    """读取模板 JSON 并检查顶层结构。"""
    if not template_json.exists():
        raise FileNotFoundError(f"模板 JSON 不存在: {template_json}")
    with template_json.open("r", encoding="utf-8") as handle:
        template = json.load(handle)
    expected_keys = ["haplogroups", "input_corrections", "metadata", "schemes"]
    if list(template.keys()) != expected_keys:
        raise PhylotreeConvertError(f"模板顶层字段不是预期结构: {list(template.keys())}")
    return template


def read_correction_sheet(input_excel: Path, sheet_name: str | None) -> dict[str, str]:
    """读取可选的旧节点名到新节点名校正表。"""
    if not sheet_name:
        return {}
    try:
        correction_df = pd.read_excel(input_excel, sheet_name=sheet_name, dtype=str)
    except ValueError:
        LOG.warning("未找到校正 sheet，跳过: %s", sheet_name)
        return {}

    required_columns = {"Parent", "Parent_Revise"}
    if not required_columns.issubset(correction_df.columns):
        LOG.warning("校正 sheet 缺少字段 %s，跳过: %s", sorted(required_columns), sheet_name)
        return {}

    corrections: dict[str, str] = {}
    for _, row in correction_df.iterrows():
        old_name = normalize_text(row.get("Parent"))
        new_name = normalize_text(row.get("Parent_Revise"))
        if old_name and new_name:
            corrections[old_name] = new_name
    return corrections


def build_output(
    input_excel: Path,
    template_json: Path,
    sheet_name: str,
    correction_sheet: str | None,
) -> dict[str, Any]:
    """生成与模板字段一致的输出对象。"""
    if not input_excel.exists():
        raise FileNotFoundError(f"输入 Excel 不存在: {input_excel}")

    template = load_template(template_json)
    LOG.info("读取 Excel sheet: %s", sheet_name)
    df = pd.read_excel(input_excel, sheet_name=sheet_name, dtype=str)
    haplogroups = build_haplogroups(df)

    input_corrections = dict(template["input_corrections"])
    input_corrections.update(read_correction_sheet(input_excel, correction_sheet))

    metadata = dict(template["metadata"])
    metadata["node_count"] = len(haplogroups)
    metadata["phylotree_table"] = f"{input_excel.name}::{sheet_name}"

    return {
        "haplogroups": haplogroups,
        "input_corrections": input_corrections,
        "metadata": metadata,
        "schemes": template["schemes"],
    }


def write_tsv(haplogroups: dict[str, dict[str, Any]], output_tsv: Path) -> None:
    """将 haplogroups 写为 TSV 文件（Haplogroup/Level/Parent/Mutations）。"""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w", encoding="utf-8", newline="") as handle:
        handle.write("Haplogroup\tLevel\tParent\tMutations\n")
        for name, record in haplogroups.items():
            handle.write(f"{name}\t{record['level']}\t{record['parent']}\t{record['mutations']}\n")
    LOG.info("TSV 已写出: %s (%d 行)", output_tsv, len(haplogroups))


def run(
    input_excel: str,
    template_json: str,
    output_json: str,
    sheet_name: str = "Phylotree build 18",
    correction_sheet: str | None = "17→18部分节点校正",
    indent: int = 2,
    output_tsv: str | None = None,
) -> dict[str, Any]:
    """执行 Excel 到 JSON 转换，并返回运行摘要。"""
    input_path = Path(input_excel)
    template_path = Path(template_json)
    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    output = build_output(input_path, template_path, sheet_name, correction_sheet)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, indent=indent)
        handle.write("\n")

    if output_tsv:
        write_tsv(output["haplogroups"], Path(output_tsv))

    summary = {
        "output_json": str(output_path),
        "node_count": len(output["haplogroups"]),
        "input_correction_count": len(output["input_corrections"]),
    }
    LOG.info("转换完成: %s", summary)
    return summary


def build_parser() -> argparse.ArgumentParser:
    """构建命令行参数解析器。"""
    parser = argparse.ArgumentParser(description="将 Phylotree build 18 Excel 转换为 JSON")
    parser.add_argument("--input-excel", required=True, help="输入 Excel 文件")
    parser.add_argument("--template-json", required=True, help="结构模板 JSON")
    parser.add_argument("--output-json", required=True, help="输出 JSON 文件")
    parser.add_argument("--sheet-name", default="Phylotree build 18", help="系统发育树 sheet 名称")
    parser.add_argument("--correction-sheet", default="17→18部分节点校正", help="节点名校正 sheet 名称")
    parser.add_argument("--indent", type=int, default=2, help="JSON 缩进空格数")
    parser.add_argument("--output-tsv", default=None, help="输出 TSV 文件（可选）")
    parser.add_argument("--log-level", default="INFO", help="日志级别")
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI 入口。"""
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    try:
        run(
            input_excel=args.input_excel,
            template_json=args.template_json,
            output_json=args.output_json,
            sheet_name=args.sheet_name,
            correction_sheet=args.correction_sheet,
            indent=args.indent,
            output_tsv=args.output_tsv,
        )
    except Exception as exc:
        LOG.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
