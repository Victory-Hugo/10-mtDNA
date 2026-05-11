#!/usr/bin/env python3
"""将已对齐 FASTA 和 TSV metadata 转换为 McAN mutation/meta/siteMask 输入。"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

log = logging.getLogger(__name__)


class McANConvertError(ValueError):
    """McAN 输入转换失败。"""


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    current_id = None
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(chunks).upper()))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                if current_id is None:
                    raise McANConvertError(f"FASTA 在第一个 header 前出现序列内容: {path}")
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, "".join(chunks).upper()))
    if not records:
        raise McANConvertError(f"FASTA 中没有读取到序列: {path}")
    return records


def read_metadata(path: Path, sample_id_column: str) -> Tuple[List[str], Dict[str, Dict[str, str]], List[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise McANConvertError(f"metadata 表头为空: {path}")
        if sample_id_column not in reader.fieldnames:
            raise McANConvertError(f"metadata 缺少样本列: {sample_id_column}")
        meta_by_id: Dict[str, Dict[str, str]] = {}
        sample_order: List[str] = []
        for row in reader:
            sample_id = (row.get(sample_id_column) or "").strip()
            if sample_id:
                meta_by_id[sample_id] = row
                sample_order.append(sample_id)
    return list(reader.fieldnames), meta_by_id, sample_order


def clean_meta_value(value: str | None, default: str) -> str:
    value = (value or "").strip()
    return value if value else default


def mutation_string(reference: str, query: str) -> Tuple[str, int]:
    mutations: List[str] = []
    skipped = 0
    for index, (ref_base, query_base) in enumerate(zip(reference, query), start=1):
        if ref_base == query_base:
            continue
        if ref_base not in "ACGT" or query_base not in "ACGT":
            skipped += 1
            continue
        mutations.append(f"{index}(SNP:{ref_base}->{query_base})")
    return ";".join(mutations), skipped


def run(
    fasta: str,
    metadata: str,
    sample_id_column: str,
    reference_id: str,
    mutation_output: str,
    meta_output: str,
    sitemask_output: str,
    report_output: str,
    country_column: str,
    state_column: str,
    city_column: str,
    default_date: str,
    unknown_label: str = "Unknown",
) -> int:
    fasta_path = Path(fasta)
    metadata_path = Path(metadata)
    for path in [fasta_path, metadata_path]:
        if not path.is_file():
            raise McANConvertError(f"输入文件不存在: {path}")

    records = read_fasta(fasta_path)
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise McANConvertError("McAN 转换要求已对齐且等长的 FASTA")
    seq_len = len(records[0][1])
    fasta_by_id = dict(records)
    if reference_id not in fasta_by_id:
        raise McANConvertError(f"reference_id 不存在于 FASTA: {reference_id}")

    fieldnames, meta_by_id, sample_order = read_metadata(metadata_path, sample_id_column)
    needed_columns = [country_column, state_column, city_column]
    missing_columns = [column for column in needed_columns if column not in fieldnames]
    if missing_columns:
        raise McANConvertError(f"metadata 缺少 McAN meta 来源列: {', '.join(missing_columns)}")
    missing_fasta = [sample_id for sample_id in sample_order if sample_id not in fasta_by_id]
    if missing_fasta:
        raise McANConvertError(f"metadata 中有样本缺失于 FASTA: {', '.join(missing_fasta[:20])}")

    ordered_ids = [reference_id] + [sample_id for sample_id in sample_order if sample_id != reference_id]
    reference_seq = fasta_by_id[reference_id]
    total_skipped = 0

    mutation_path = Path(mutation_output)
    meta_path = Path(meta_output)
    sitemask_path = Path(sitemask_output)
    report_path = Path(report_output)
    for path in [mutation_path, meta_path, sitemask_path, report_path]:
        path.parent.mkdir(parents=True, exist_ok=True)

    with mutation_path.open("w", encoding="utf-8", newline="\n") as mut_handle:
        for sample_id in ordered_ids:
            if sample_id == reference_id:
                mut_handle.write(f"{sample_id}\t{sample_id}\n")
                continue
            mutations, skipped = mutation_string(reference_seq, fasta_by_id[sample_id])
            total_skipped += skipped
            if mutations:
                mut_handle.write(f"{sample_id}\t{sample_id}\t{mutations}\n")
            else:
                mut_handle.write(f"{sample_id}\t{sample_id}\n")

    with meta_path.open("w", encoding="utf-8", newline="\n") as meta_handle:
        for sample_id in ordered_ids:
            row = meta_by_id.get(sample_id, {})
            country = clean_meta_value(row.get(country_column), unknown_label)
            state = clean_meta_value(row.get(state_column), unknown_label)
            city = clean_meta_value(row.get(city_column), unknown_label)
            meta_handle.write(f"{sample_id}\t{sample_id}\t{default_date}\t{country}\t{state}\t{city}\n")

    with sitemask_path.open("w", encoding="utf-8", newline="\n") as mask_handle:
        for pos in range(1, seq_len + 1):
            mask_handle.write(f"{pos}\n")

    with report_path.open("w", encoding="utf-8", newline="\n") as report_handle:
        report_handle.write("metric\tvalue\n")
        report_handle.write(f"reference_id\t{reference_id}\n")
        report_handle.write(f"sequence_length\t{seq_len}\n")
        report_handle.write(f"samples_written\t{len(ordered_ids)}\n")
        report_handle.write(f"skipped_non_acgt_differences\t{total_skipped}\n")

    log.info("McAN 输入转换完成: %s 条样本", len(ordered_ids))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="已过滤排序的 FASTA")
    parser.add_argument("--metadata", required=True, help="metadata TSV")
    parser.add_argument("--sample-id-column", required=True, help="样本 ID 列名")
    parser.add_argument("--reference-id", required=True, help="McAN 参考序列 ID")
    parser.add_argument("--mutation-output", required=True, help="输出 McAN mutation 文件")
    parser.add_argument("--meta-output", required=True, help="输出 McAN meta 文件")
    parser.add_argument("--sitemask-output", required=True, help="输出 McAN siteMask 文件")
    parser.add_argument("--report-output", required=True, help="输出转换报告")
    parser.add_argument("--country-column", required=True, help="Country 来源列")
    parser.add_argument("--state-column", required=True, help="State 来源列")
    parser.add_argument("--city-column", required=True, help="City 来源列")
    parser.add_argument("--default-date", required=True, help="McAN 占位采样日期")
    parser.add_argument("--unknown-label", default="Unknown", help="缺失值标签")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(**vars(args))
    except McANConvertError as exc:
        parser.error(str(exc))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
