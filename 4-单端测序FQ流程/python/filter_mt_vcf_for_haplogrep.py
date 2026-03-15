#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def run(
    input: str,
    output: str,
    stats: str,
    min_dp: int,
    min_vaf: float,
    pass_only: bool,
) -> dict[str, int | float | str]:
    in_path = Path(input)
    out_path = Path(output)
    stats_path = Path(stats)
    if not in_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {in_path}")

    total_records = 0
    kept_records = 0
    filtered_by_status = 0
    filtered_by_dp = 0
    filtered_by_vaf = 0

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(in_path, "rt") as reader, open(out_path, "w", encoding="utf-8") as writer:
        for line in reader:
            if line.startswith("#"):
                writer.write(line)
                continue

            total_records += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            filt = fields[6]
            if pass_only and filt not in {"PASS", "."}:
                filtered_by_status += 1
                continue

            fmt_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            fmt = dict(zip(fmt_keys, sample_values))

            try:
                dp = int(fmt.get("DP", "0") or "0")
            except ValueError:
                dp = 0
            if dp < min_dp:
                filtered_by_dp += 1
                continue

            ad_raw = fmt.get("AD", "")
            if not ad_raw:
                filtered_by_vaf += 1
                continue

            try:
                ad_values = [int(x) for x in ad_raw.split(",") if x != "."]
            except ValueError:
                filtered_by_vaf += 1
                continue

            alt_depth = max(ad_values[1:], default=0)
            vaf = alt_depth / dp if dp else 0.0
            if vaf < min_vaf:
                filtered_by_vaf += 1
                continue

            kept_records += 1
            writer.write(line)

    results = {
        "input_vcf": str(in_path),
        "output_vcf": str(out_path),
        "total_records": total_records,
        "kept_records": kept_records,
        "filtered_by_status": filtered_by_status,
        "filtered_by_dp": filtered_by_dp,
        "filtered_by_vaf": filtered_by_vaf,
        "min_dp": min_dp,
        "min_vaf": min_vaf,
        "pass_only": int(pass_only),
    }
    with open(stats_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in results.items():
            writer.writerow([key, value])
    log.info("Filtered VCF written to %s", out_path)
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Filter a single-sample mtDNA VCF for Haplogrep.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--min-dp", type=int, required=True)
    parser.add_argument("--min-vaf", type=float, required=True)
    parser.add_argument("--pass-only", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
