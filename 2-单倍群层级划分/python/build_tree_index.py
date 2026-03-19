#!/usr/bin/env python3
"""Build a JSON index for mtDNA haplogroup annotation."""

from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
from typing import Dict, Iterable, List

log = logging.getLogger(__name__)


def normalize_key(name: str) -> str:
    return name.strip().lower().replace(" ", "").replace("_", "")


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)s] %(message)s",
    )


def read_mapping_table(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if len(row) < 2:
                continue
            source = row[0].strip()
            target = row[1].strip()
            if not source or not target:
                continue
            if source.lower() == "origin" and target.lower() == "revised":
                continue
            mapping[source] = target
    return mapping


def read_target_list(path: Path) -> List[str]:
    targets: List[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for raw_line in handle:
            value = raw_line.strip()
            if value:
                targets.append(value)
    return targets


def read_phylotree_table(path: Path) -> Dict[str, Dict[str, object]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in PhyloTree table: {path}")
        field_map = {normalize_key(name): name for name in reader.fieldnames}
        required_keys = {
            "haplogroup": "Haplogroup",
            "level": "Level",
            "parent": "Parent",
            "mutations": "Mutations",
        }
        missing = [alias for alias in required_keys if alias not in field_map]
        if missing:
            raise ValueError(f"Missing columns in {path}: {', '.join(missing)}")

        nodes: Dict[str, Dict[str, object]] = {}
        for row in reader:
            haplogroup = row[field_map["haplogroup"]].strip()
            if not haplogroup:
                continue
            level_text = row[field_map["level"]].strip()
            parent = row[field_map["parent"]].strip()
            mutations = row[field_map["mutations"]].strip()
            nodes[haplogroup] = {
                "level": int(level_text),
                "parent": parent,
                "mutations": mutations,
            }
    return nodes


def build_lineage_map(nodes: Dict[str, Dict[str, object]]) -> Dict[str, List[str]]:
    lineage_cache: Dict[str, List[str]] = {}

    def lineage_for(name: str) -> List[str]:
        if name in lineage_cache:
            return lineage_cache[name]

        node = nodes.get(name)
        if node is None:
            raise KeyError(f"Unknown haplogroup in lineage lookup: {name}")

        parent = str(node["parent"])
        if not parent:
            lineage = [name]
        else:
            if parent not in nodes:
                raise KeyError(f"Parent haplogroup missing from tree: {parent}")
            lineage = lineage_for(parent) + [name]
        lineage_cache[name] = lineage
        return lineage

    for haplogroup in nodes:
        lineage_for(haplogroup)

    return lineage_cache


def run(
    phylotree_table: str,
    input_correction: str,
    llt_targets: str,
    llt_correction: str,
    yuchunli_targets: str,
    yuchunli_correction: str,
    output: str,
    log_level: str = "INFO",
) -> int:
    configure_logging(log_level)

    phylotree_path = Path(phylotree_table)
    output_path = Path(output)

    log.info("Loading PhyloTree table: %s", phylotree_path)
    nodes = read_phylotree_table(phylotree_path)
    log.info("Loaded %d PhyloTree haplogroups", len(nodes))

    lineages = build_lineage_map(nodes)

    data = {
        "metadata": {
            "phylotree_table": str(phylotree_path),
            "node_count": len(nodes),
        },
        "input_corrections": read_mapping_table(Path(input_correction)),
        "schemes": {
            "LLT": {
                "targets": read_target_list(Path(llt_targets)),
                "corrections": read_mapping_table(Path(llt_correction)),
            },
            "YuChunLi": {
                "targets": read_target_list(Path(yuchunli_targets)),
                "corrections": read_mapping_table(Path(yuchunli_correction)),
            },
        },
        "haplogroups": {
            name: {
                **details,
                "lineage": lineages[name],
            }
            for name, details in nodes.items()
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2, sort_keys=True)

    log.info("Tree index written to %s", output_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phylotree-table", required=True, help="PhyloTree TSV file.")
    parser.add_argument("--input-correction", required=True, help="Input haplogroup correction TSV.")
    parser.add_argument("--llt-targets", required=True, help="LLT target haplogroup list.")
    parser.add_argument("--llt-correction", required=True, help="LLT correction TSV.")
    parser.add_argument("--yuchunli-targets", required=True, help="YuChunLi target haplogroup list.")
    parser.add_argument("--yuchunli-correction", required=True, help="YuChunLi correction TSV.")
    parser.add_argument("--output", required=True, help="Output JSON index path.")
    parser.add_argument("--log-level", default="INFO", help="Logging level.")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
