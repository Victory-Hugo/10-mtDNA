#!/usr/bin/env python3
"""解析 PhyloTree Build 17 HTM 文件，输出层级关系和定义突变。"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from typing import Iterable, List

from bs4 import BeautifulSoup

log = logging.getLogger(__name__)

HAPLO_CLASSES = {
    "xl11317826",
    "xl11417826",
    "xl11517826",
    "xl11617826",
    "xl11717826",
    "xl11817826",
    "xl11917826",
    "xl13417826",
}
ACCESSION_CLASS = "xl13117826"
MUTATION_CLASS_PREFIXES = ("xl93", "xl97")


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)s] %(message)s",
    )


def clean_mutation_text(text: str) -> str:
    text = text.replace("\ufeff", "").replace("\xa0", "").replace("\u00a0", "")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def extract_mutations_from_row(cells_after_haplo) -> str:
    parts = []
    for td in cells_after_haplo:
        cls_list = td.get("class", [])
        cls_str = " ".join(cls_list) if isinstance(cls_list, list) else cls_list
        if ACCESSION_CLASS in cls_str:
            break
        if any(cls_str.startswith(prefix) or prefix in cls_str for prefix in MUTATION_CLASS_PREFIXES):
            text = clean_mutation_text(td.get_text(" "))
            if text:
                parts.append(text)
    return " ".join(parts)


def is_haplo_cell(td) -> bool:
    if td.get("rowspan") != "2":
        return False
    cls_list = td.get("class", [])
    cls_set = set(cls_list) if isinstance(cls_list, list) else {cls_list}
    return bool(cls_set & HAPLO_CLASSES)


def parse_phylotree(filepath: str) -> List[dict]:
    path = Path(filepath)
    log.info("读取文件 %s", path)
    content = path.read_text(encoding="windows-1252", errors="replace")

    soup = BeautifulSoup(content, "lxml")
    table = soup.find("table")
    if not table:
        raise ValueError("找不到 <table>")

    rows = table.find_all("tr")
    log.info("HTML 行数: %d", len(rows))

    haplogroups = []
    depth_stack = {}

    for tr in rows:
        cells = list(tr.find_all("td"))
        if not cells:
            continue

        col = 0
        haplo_name = None
        haplo_col = None
        cells_after = []

        for idx, td in enumerate(cells):
            colspan = int(td.get("colspan", 1))
            if haplo_name is None and is_haplo_cell(td):
                text = clean_mutation_text(td.get_text(" "))
                if text:
                    haplo_name = text
                    haplo_col = col
                    cells_after = cells[idx + 1 :]
            col += colspan

        if not haplo_name:
            continue
        if "mt-MRCA" in haplo_name or "MRCA" in haplo_name:
            depth_stack[haplo_col] = haplo_name
            continue

        parent = ""
        for depth in range(haplo_col - 1, -1, -1):
            if depth in depth_stack:
                parent = depth_stack[depth]
                break

        depth_stack[haplo_col] = haplo_name
        for depth in list(depth_stack.keys()):
            if depth > haplo_col:
                del depth_stack[depth]

        haplogroups.append(
            {
                "haplogroup": haplo_name,
                "depth": haplo_col,
                "parent": parent,
                "mutations": extract_mutations_from_row(cells_after),
            }
        )

    return haplogroups


def write_output(rows: List[dict], output: str | None) -> None:
    fieldnames = ["haplogroup", "depth", "parent", "mutations"]
    if output:
        out_path = Path(output)
        with out_path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        log.info("已写入：%s", out_path)
        return

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


def run(input: str, output: str | None = None, log_level: str = "INFO") -> int:
    configure_logging(log_level)
    rows = parse_phylotree(input)
    log.info("共提取 %d 个单倍群节点", len(rows))
    write_output(rows, output)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="PhyloTree HTM 文件")
    parser.add_argument("--output", help="输出 CSV 路径，不给则写到 stdout")
    parser.add_argument("--log-level", default="INFO", help="日志级别")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
