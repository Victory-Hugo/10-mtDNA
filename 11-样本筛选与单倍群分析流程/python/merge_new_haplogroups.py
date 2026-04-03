#!/usr/bin/env python3
"""把 txt 文件中含下划线的新单倍群合并到 phylotree_build17 CSV。"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

log = logging.getLogger(__name__)


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)s] %(message)s",
    )


def read_txt_nodes(txt_path: Path) -> Tuple[List[Tuple[int, str]], Dict[str, str]]:
    txt_nodes: List[Tuple[int, str]] = []
    with txt_path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            name = line.lstrip("\t")
            indent = len(line) - len(name)
            name = name.strip()
            if name:
                txt_nodes.append((indent, name))

    parent_of: Dict[str, str] = {}
    stack: List[Tuple[int, str]] = []
    for indent, name in txt_nodes:
        while stack and stack[-1][0] >= indent:
            stack.pop()
        parent_of[name] = stack[-1][1] if stack else ""
        stack.append((indent, name))
    return txt_nodes, parent_of


def merge_rows(txt_nodes: List[Tuple[int, str]], parent_of: Dict[str, str], csv_path: Path) -> List[dict]:
    existing_rows = []
    csv_level = {}
    with csv_path.open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            existing_rows.append(row)
            csv_level[row["haplogroup"]] = int(row["level"])

    csv_names = set(csv_level)
    memo: Dict[str, int] = {}

    def compute_level(name: str) -> int:
        if name in memo:
            return memo[name]
        if name in csv_level:
            memo[name] = csv_level[name]
            return csv_level[name]
        parent = parent_of.get(name, "")
        memo[name] = 0 if not parent else compute_level(parent) + 1
        return memo[name]

    new_haplos = [name for _, name in txt_nodes if "_" in name and name not in csv_names]
    log.info("含下划线新单倍群：%d", len(new_haplos))

    new_row_of = {}
    for name in new_haplos:
        parts = name.split("_")
        new_row_of[name] = {
            "haplogroup": name,
            "level": compute_level(name),
            "parent": parent_of.get(name, ""),
            "mutations": " ".join(parts[1:]),
        }

    csv_row_of = {row["haplogroup"]: row for row in existing_rows}
    output_rows = []
    seen = set()

    for _, name in txt_nodes:
        if name in seen:
            continue
        seen.add(name)
        if name in csv_row_of:
            output_rows.append(csv_row_of[name])
        elif name in new_row_of:
            output_rows.append(new_row_of[name])

    covered = {row["haplogroup"] for row in output_rows}
    not_covered = [name for name in csv_names if name not in covered]
    if not_covered:
        log.warning("%d 个 CSV 节点未出现在 txt 树中，追加到末尾", len(not_covered))
        for name in not_covered:
            output_rows.append(csv_row_of[name])

    return output_rows


def run(input_txt: str, input_csv: str, output_csv: str, log_level: str = "INFO") -> int:
    configure_logging(log_level)
    txt_path = Path(input_txt)
    csv_path = Path(input_csv)
    output_path = Path(output_csv)

    txt_nodes, parent_of = read_txt_nodes(txt_path)
    log.info("txt 节点总数：%d", len(txt_nodes))
    output_rows = merge_rows(txt_nodes, parent_of, csv_path)

    with output_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["haplogroup", "level", "parent", "mutations"], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)
    log.info("完成，已写入 %s", output_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-txt", required=True, help="树结构 txt 文件")
    parser.add_argument("--input-csv", required=True, help="现有 phylotree CSV/TSV 文件")
    parser.add_argument("--output-csv", required=True, help="合并后的输出 CSV")
    parser.add_argument("--log-level", default="INFO", help="日志级别")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
