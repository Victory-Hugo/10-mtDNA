#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def count_vcf_records(vcf_path: Path) -> int:
    count = 0
    with open(vcf_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                count += 1
    return count


def parse_haplogrep_top_hit(path: Path) -> dict[str, str]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No Haplogrep rows found: {path}")
    return rows[0]


def parse_filter_stats(path: Path) -> dict[str, str]:
    stats: dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            stats[row["metric"]] = row["value"]
    return stats


def parse_coverage(path: Path) -> dict[str, str]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No coverage rows found: {path}")
    return rows[0]


def run(
    sample_id: str,
    raw_vcf: str,
    strict_vcf: str,
    raw_haplogrep: str,
    strict_haplogrep: str,
    strict_filter_stats: str,
    coverage_tsv: str,
    output_md: str,
) -> str:
    raw_vcf_path = Path(raw_vcf)
    strict_vcf_path = Path(strict_vcf)
    raw_haplogrep_path = Path(raw_haplogrep)
    strict_haplogrep_path = Path(strict_haplogrep)
    strict_filter_stats_path = Path(strict_filter_stats)
    coverage_path = Path(coverage_tsv)

    raw_count = count_vcf_records(raw_vcf_path)
    strict_count = count_vcf_records(strict_vcf_path)
    raw_hit = parse_haplogrep_top_hit(raw_haplogrep_path)
    strict_hit = parse_haplogrep_top_hit(strict_haplogrep_path)
    filter_stats = parse_filter_stats(strict_filter_stats_path)
    coverage = parse_coverage(coverage_path)

    lines = [
        f"# {sample_id} summary",
        "",
        "## Coverage",
        f"- chrM coverage breadth: {coverage['coverage']}%",
        f"- chrM mean depth: {coverage['meandepth']}x",
        "",
        "## VCF",
        f"- Raw VCF records: {raw_count}",
        f"- Strict VCF records: {strict_count}",
        f"- Strict filter threshold: DP>={filter_stats['min_dp']}, VAF>={filter_stats['min_vaf']}",
        f"- Variants removed by strict filter: {int(filter_stats['total_records']) - int(filter_stats['kept_records'])}",
        "",
        "## Haplogrep3",
        f"- Raw top hit: {raw_hit['Haplogroup']} (rank {raw_hit['Rank']}, quality {raw_hit['Quality']})",
        f"- Strict top hit: {strict_hit['Haplogroup']} (rank {strict_hit['Rank']}, quality {strict_hit['Quality']})",
        "",
        "## Interpretation",
    ]
    if raw_hit["Haplogroup"] == strict_hit["Haplogroup"]:
        lines.append(
            f"- Raw and strict VCFs agree on haplogroup `{raw_hit['Haplogroup']}`, so the major lineage signal is stable."
        )
    else:
        lines.append(
            f"- Raw and strict VCFs disagree (`{raw_hit['Haplogroup']}` vs `{strict_hit['Haplogroup']}`), so this sample should be reviewed manually."
        )
    lines.append(
        f"- Raw quality is {raw_hit['Quality']}; strict quality is {strict_hit['Quality']}."
    )
    summary = "\n".join(lines) + "\n"

    output_path = Path(output_md)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(summary, encoding="utf-8")
    log.info("Summary written to %s", output_path)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Summarize mtDNA VCF and Haplogrep results.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--raw-vcf", required=True)
    parser.add_argument("--strict-vcf", required=True)
    parser.add_argument("--raw-haplogrep", required=True)
    parser.add_argument("--strict-haplogrep", required=True)
    parser.add_argument("--strict-filter-stats", required=True)
    parser.add_argument("--coverage-tsv", required=True)
    parser.add_argument("--output-md", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
