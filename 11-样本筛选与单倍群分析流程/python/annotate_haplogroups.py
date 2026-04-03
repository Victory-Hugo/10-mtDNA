#!/usr/bin/env python3
"""Annotate mtDNA haplogroups against a PhyloTree-derived JSON index."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

log = logging.getLogger(__name__)

OUTPUT_COLUMNS = [
    "ID",
    "Haplogroup_Original",
    "Haplogroup_Standardized",
    "Haplogroup_PhyloTree",
    "Haplogroup_LLT",
    "Haplogroup_YuChunLi",
    "Resolution_Status",
    "Phylotree_Level",
    "Phylotree_Parent",
    "Phylotree_Mutations",
]


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)s] %(message)s",
    )


def load_tree_index(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def detect_dialect(sample: str) -> csv.Dialect:
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        class DefaultTabDialect(csv.excel_tab):
            delimiter = "\t"

        return DefaultTabDialect()


def looks_like_header(first_row: Sequence[str], sniffer_guess: bool) -> bool:
    if len(first_row) < 2:
        return sniffer_guess

    first = first_row[0].strip().lower()
    second = first_row[1].strip().lower()
    header_tokens = {"id", "sample", "sampleid", "sample_id"}
    if first in header_tokens:
        return True
    if "hap" in second:
        return True
    return sniffer_guess


def read_input_rows(path: Path) -> List[Tuple[str, str]]:
    text = path.read_text(encoding="utf-8-sig")
    if not text.strip():
        return []

    sample = text[:4096]
    dialect = detect_dialect(sample)
    sniffer = csv.Sniffer()
    try:
        header_guess = sniffer.has_header(sample)
    except csv.Error:
        header_guess = False

    rows: List[List[str]] = []
    for row in csv.reader(text.splitlines(), dialect):
        if not row:
            continue
        trimmed = [value.strip() for value in row]
        if not any(trimmed):
            continue
        rows.append(trimmed)

    if not rows:
        return []

    start_index = 1 if looks_like_header(rows[0], header_guess) else 0
    result: List[Tuple[str, str]] = []
    for row in rows[start_index:]:
        if len(row) < 2:
            raise ValueError(f"Input row has fewer than 2 columns: {row}")
        result.append((row[0], row[1]))
    return result


def infer_from_prefix(candidate: str, nodes: Dict[str, Dict[str, object]]) -> str:
    value = candidate.strip()
    while value:
        if value in nodes:
            return value
        value = value[:-1]
    return ""


def resolve_scheme(lineage: Sequence[str], scheme: Dict[str, object]) -> str:
    targets = set(scheme["targets"])
    corrections = dict(scheme["corrections"])
    for haplogroup in reversed(lineage):
        if haplogroup in targets:
            return corrections.get(haplogroup, haplogroup)
    return ""


def resolve_haplogroup(
    original_haplogroup: str,
    input_corrections: Dict[str, str],
    nodes: Dict[str, Dict[str, object]],
) -> Tuple[str, str, str]:
    standardized = input_corrections.get(original_haplogroup, original_haplogroup)
    if standardized in nodes:
        if standardized == original_haplogroup:
            return standardized, standardized, "exact"
        return standardized, standardized, "corrected"

    inferred = infer_from_prefix(standardized, nodes)
    if inferred:
        if standardized == original_haplogroup:
            return standardized, inferred, "inferred_ancestor"
        return standardized, inferred, "corrected_inferred_ancestor"

    return standardized, "", "unresolved"


def build_output_rows(
    input_rows: Sequence[Tuple[str, str]],
    tree_index: Dict[str, object],
) -> List[Dict[str, object]]:
    nodes = dict(tree_index["haplogroups"])
    input_corrections = dict(tree_index["input_corrections"])
    schemes = dict(tree_index["schemes"])

    output_rows: List[Dict[str, object]] = []
    unresolved_count = 0

    for sample_id, original_haplogroup in input_rows:
        standardized, resolved, status = resolve_haplogroup(
            original_haplogroup,
            input_corrections,
            nodes,
        )
        node = nodes.get(resolved, {})
        lineage = list(node.get("lineage", []))
        llt_label = resolve_scheme(lineage, schemes["LLT"]) if lineage else ""
        yuchunli_label = resolve_scheme(lineage, schemes["YuChunLi"]) if lineage else ""

        if status == "unresolved":
            unresolved_count += 1

        output_rows.append(
            {
                "ID": sample_id,
                "Haplogroup_Original": original_haplogroup,
                "Haplogroup_Standardized": standardized,
                "Haplogroup_PhyloTree": resolved,
                "Haplogroup_LLT": llt_label,
                "Haplogroup_YuChunLi": yuchunli_label,
                "Resolution_Status": status,
                "Phylotree_Level": node.get("level", ""),
                "Phylotree_Parent": node.get("parent", ""),
                "Phylotree_Mutations": node.get("mutations", ""),
            }
        )

    log.info("Annotated %d rows", len(output_rows))
    if unresolved_count:
        log.warning("Unresolved haplogroups: %d", unresolved_count)
    return output_rows


def write_output(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run(input: str, tree_index: str, output: str, log_level: str = "INFO") -> int:
    configure_logging(log_level)

    input_path = Path(input)
    tree_index_path = Path(tree_index)
    output_path = Path(output)

    log.info("Reading input haplogroups: %s", input_path)
    input_rows = read_input_rows(input_path)

    log.info("Loading tree index: %s", tree_index_path)
    tree_data = load_tree_index(tree_index_path)

    output_rows = build_output_rows(input_rows, tree_data)
    write_output(output_path, output_rows)
    log.info("Annotation written to %s", output_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input TSV or CSV file.")
    parser.add_argument("--tree-index", required=True, help="JSON tree index file.")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
