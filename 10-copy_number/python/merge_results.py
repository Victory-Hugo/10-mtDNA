#!/usr/bin/env python3
import argparse
import csv
import os
from typing import List


def read_sample_rows(path: str) -> List[List[str]]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty file: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No data rows in: {path}")
    return rows


def run(output_dir: str, output_path: str) -> None:
    files = [
        os.path.join(output_dir, name)
        for name in os.listdir(output_dir)
        if name.endswith(".tsv") and os.path.join(output_dir, name) != output_path
    ]
    if not files:
        raise SystemExit("No per-sample TSV files found")

    all_rows: List[List[str]] = []
    for path in sorted(files):
        all_rows.extend(read_sample_rows(path))

    all_rows.sort(key=lambda row: row[0])

    tmp_path = f"{output_path}.tmp"
    with open(tmp_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "mean_mtDNA_depth",
                "mean_autosomal_depth",
                "mtDNA_copy_number",
            ]
        )
        writer.writerows(all_rows)
    os.replace(tmp_path, output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge per-sample mtDNA copy number TSVs.")
    parser.add_argument("--output-dir", required=True, help="Directory with per-sample TSVs")
    parser.add_argument("--output", required=True, help="Merged TSV path")
    args = parser.parse_args()
    run(args.output_dir, args.output)


if __name__ == "__main__":
    main()
