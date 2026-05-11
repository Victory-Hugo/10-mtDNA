#!/usr/bin/env python3
"""校验 network pipeline 输入，并生成 fastHaN 使用的 PHYLIP 文件。"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

log = logging.getLogger(__name__)


class InputError(ValueError):
    """输入数据不满足 pipeline 要求。"""


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
                    raise InputError(f"FASTA 在第一个 header 前出现序列内容: {path}")
                chunks.append(line)

    if current_id is not None:
        records.append((current_id, "".join(chunks).upper()))
    if not records:
        raise InputError(f"FASTA 中没有读取到序列: {path}")
    return records


def read_metadata(path: Path, sample_id_column: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise InputError(f"metadata 表头为空: {path}")
        if sample_id_column not in reader.fieldnames:
            raise InputError(f"metadata 缺少样本列 {sample_id_column}: {path}")
        rows = [dict(row) for row in reader]

    sample_ids = [(row.get(sample_id_column) or "").strip() for row in rows]
    if any(not sample_id for sample_id in sample_ids):
        raise InputError(f"metadata 的 {sample_id_column} 列存在空值")
    return list(reader.fieldnames), rows


def duplicate_items(items: Sequence[str]) -> List[str]:
    seen = set()
    dup = []
    for item in items:
        if item in seen and item not in dup:
            dup.append(item)
        seen.add(item)
    return dup


def write_fasta(records: Iterable[Tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for sample_id, sequence in records:
            handle.write(f">{sample_id}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")


def write_phylip(records: Sequence[Tuple[str, str]], path: Path) -> None:
    if not records:
        raise InputError("没有可写入 PHYLIP 的序列")
    seq_len = len(records[0][1])
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{len(records)} {seq_len}\n")
        for sample_id, sequence in records:
            handle.write(f"{sample_id} {sequence}\n")


def write_qc_report(
    path: Path,
    fasta_records: Sequence[Tuple[str, str]],
    metadata_rows: Sequence[Dict[str, str]],
    sample_id_column: str,
    ordered_records: Sequence[Tuple[str, str]],
) -> None:
    fasta_ids = [sample_id for sample_id, _ in fasta_records]
    meta_ids = [(row.get(sample_id_column) or "").strip() for row in metadata_rows]
    fasta_id_set = set(fasta_ids)
    meta_id_set = set(meta_ids)
    non_acgt = sum(sum(1 for base in sequence if base not in "ACGT") for _, sequence in ordered_records)
    unique_lengths = sorted({len(sequence) for _, sequence in ordered_records})

    rows = [
        ("fasta_records", len(fasta_records)),
        ("metadata_records", len(metadata_rows)),
        ("matched_records", len(ordered_records)),
        ("sequence_lengths", ",".join(str(x) for x in unique_lengths)),
        ("missing_metadata_ids_in_fasta", len(meta_id_set - fasta_id_set)),
        ("extra_fasta_ids_not_in_metadata", len(fasta_id_set - meta_id_set)),
        ("duplicate_fasta_ids", len(duplicate_items(fasta_ids))),
        ("duplicate_metadata_ids", len(duplicate_items(meta_ids))),
        ("non_acgt_bases_in_matched_sequences", non_acgt),
    ]
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("metric\tvalue\n")
        for key, value in rows:
            handle.write(f"{key}\t{value}\n")


def run(
    fasta: str,
    metadata: str,
    sample_id_column: str,
    filtered_fasta: str,
    phylip: str,
    qc_report: str,
) -> int:
    fasta_path = Path(fasta)
    metadata_path = Path(metadata)
    filtered_fasta_path = Path(filtered_fasta)
    phylip_path = Path(phylip)
    qc_report_path = Path(qc_report)

    for path in [fasta_path, metadata_path]:
        if not path.is_file():
            raise InputError(f"输入文件不存在: {path}")

    fasta_records = read_fasta(fasta_path)
    _fieldnames, metadata_rows = read_metadata(metadata_path, sample_id_column)
    fasta_ids = [sample_id for sample_id, _ in fasta_records]
    meta_ids = [(row.get(sample_id_column) or "").strip() for row in metadata_rows]

    duplicate_fasta = duplicate_items(fasta_ids)
    duplicate_meta = duplicate_items(meta_ids)
    if duplicate_fasta:
        raise InputError(f"FASTA 存在重复样本 ID: {', '.join(duplicate_fasta[:10])}")
    if duplicate_meta:
        raise InputError(f"metadata 存在重复样本 ID: {', '.join(duplicate_meta[:10])}")

    fasta_by_id = dict(fasta_records)
    missing = [sample_id for sample_id in meta_ids if sample_id not in fasta_by_id]
    if missing:
        raise InputError(f"metadata 中有样本缺失于 FASTA: {', '.join(missing[:20])}")

    ordered_records = [(sample_id, fasta_by_id[sample_id]) for sample_id in meta_ids]
    lengths = {len(sequence) for _, sequence in ordered_records}
    if len(lengths) != 1:
        preview = ", ".join(str(x) for x in sorted(lengths)[:10])
        raise InputError(f"FASTA 序列长度不一致，当前 pipeline 要求已对齐输入。长度: {preview}")

    filtered_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    phylip_path.parent.mkdir(parents=True, exist_ok=True)
    qc_report_path.parent.mkdir(parents=True, exist_ok=True)

    write_fasta(ordered_records, filtered_fasta_path)
    write_phylip(ordered_records, phylip_path)
    write_qc_report(qc_report_path, fasta_records, metadata_rows, sample_id_column, ordered_records)
    log.info("输入整理完成: %s 条序列", len(ordered_records))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="输入 FASTA 文件")
    parser.add_argument("--metadata", required=True, help="输入 metadata TSV")
    parser.add_argument("--sample-id-column", required=True, help="样本 ID 列名")
    parser.add_argument("--filtered-fasta", required=True, help="输出过滤排序后的 FASTA")
    parser.add_argument("--phylip", required=True, help="输出 PHYLIP 文件")
    parser.add_argument("--qc-report", required=True, help="输出 QC 报告 TSV")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run(**vars(args))
    except InputError as exc:
        parser.error(str(exc))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
