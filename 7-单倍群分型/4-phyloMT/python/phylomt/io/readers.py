"""HSD、VCF、FASTA 输入解析。"""

from __future__ import annotations

from pathlib import Path

import pysam
from Bio import SeqIO

from ..models import SampleProfile
from ..mutations import apply_alignment_rules, clean_token, ordered_unique


def load_samples(
    input_path: Path,
    reference_sequence: str,
    rules: list[tuple[list[str], list[str]]],
    het_level: float,
    threads: int = 1,
    chip: bool = False,
) -> list[SampleProfile]:
    """按输入格式分发样本解析。"""

    suffixes = "".join(input_path.suffixes).lower()
    if suffixes.endswith(".vcf") or suffixes.endswith(".vcf.gz"):
        return read_vcf(input_path, het_level, chip=chip)
    if suffixes.endswith(".fasta") or suffixes.endswith(".fa") or suffixes.endswith(".fasta.gz") or suffixes.endswith(".fa.gz"):
        return read_fasta(input_path, reference_sequence, rules, threads=threads)
    return read_hsd(input_path)


def read_hsd(input_path: Path) -> list[SampleProfile]:
    """读取 HSD 文本（Tab 分隔）。"""

    from ..mutations import parse_range_to_positions

    samples: list[SampleProfile] = []
    with input_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            if fields[0].strip().lower() in {"id", "sampleid", "sample"}:
                continue
            sample_id = fields[0].strip()
            raw_range = fields[1].strip().strip('"').strip()
            range_text = " ".join(p.strip() for p in raw_range.split(";") if p.strip())
            covered_positions = parse_range_to_positions(range_text)
            variants = ordered_unique([
                clean_token(token) for token in fields[3:]
                if clean_token(token) and not clean_token(token).endswith("R")
            ])
            samples.append(
                SampleProfile(
                    sample_id=sample_id,
                    range_text=range_text,
                    variants=variants,
                    source_path=input_path,
                    covered_positions=covered_positions,
                )
            )
    return samples


def read_vcf(input_path: Path, het_level: float, chip: bool = False) -> list[SampleProfile]:
    """读取 VCF 文件。"""

    sample_variants: dict[str, list[str]] = {}
    covered_pos: set[int] = set()

    with pysam.VariantFile(str(input_path)) as variant_file:
        sample_names = list(variant_file.header.samples)
        for sample_name in sample_names:
            sample_variants[sample_name] = []

        for record in variant_file.fetch():
            if chip:
                covered_pos.add(record.pos)
            for sample_name in sample_names:
                call = record.samples[sample_name]
                if not _sample_has_alt(call, het_level):
                    continue
                sample_variants[sample_name].extend(vcf_record_to_tokens(record, call))

    covered_positions = covered_pos if chip else None
    return [
        SampleProfile(
            sample_id=sample_name,
            range_text="1-16569",
            variants=ordered_unique(tokens),
            source_path=input_path,
            covered_positions=covered_positions,
        )
        for sample_name, tokens in sample_variants.items()
    ]


def _sample_has_alt(call: pysam.libcbcf.VariantRecordSample, het_level: float) -> bool:
    """判断样本是否携带 ALT。"""

    gt = call.get("GT")
    if not gt:
        return False
    alt_count = sum(allele not in (None, 0) for allele in gt)
    if alt_count == 0:
        return False
    ad = call.get("AD")
    if ad and len(ad) > 1:
        depth = sum(ad)
        if depth > 0:
            alt_fraction = sum(ad[1:]) / depth
            return alt_fraction >= het_level
    return True


def vcf_record_to_tokens(record: pysam.VariantRecord, call: "pysam.libcbcf.VariantRecordSample | None" = None) -> list[str]:
    """把单条 VCF 记录转为变异标记（只处理样本 GT 中实际携带的 ALT）。"""

    tokens: list[str] = []
    ref = record.ref.upper()
    all_alts = [alt.upper() for alt in record.alts or []]
    if not all_alts:
        return tokens

    if call is not None:
        gt = call.get("GT") or ()
        alt_indices = {allele - 1 for allele in gt if allele is not None and allele > 0}
        alts = [all_alts[i] for i in sorted(alt_indices) if i < len(all_alts)]
    else:
        alts = all_alts

    for alt in alts:
        if len(ref) == 1 and len(alt) == 1:
            tokens.append(f"{record.pos}{alt}")
            continue
        if len(ref) == 1 and len(alt) > 1:
            inserted = alt[1:]
            for index, base in enumerate(inserted, start=1):
                tokens.append(f"{record.pos}.{index}{base}")
            continue
        if len(ref) > 1 and len(alt) == 1:
            for offset in range(1, len(ref)):
                tokens.append(f"{record.pos + offset}d")
            continue
        if len(ref) == len(alt):
            for offset, (ref_base, alt_base) in enumerate(zip(ref, alt)):
                if ref_base != alt_base:
                    tokens.append(f"{record.pos + offset}{alt_base}")
            continue
    return tokens


def read_fasta(
    input_path: Path,
    reference_sequence: str,
    rules: list[tuple[list[str], list[str]]],
    threads: int = 1,
) -> list[SampleProfile]:
    """读取 FASTA 并转换为变异列表（使用 halign4 批量比对）。"""

    from ..alignment import align_batch_with_halign4, normalize_circular_query, variants_from_alignment

    records = list(SeqIO.parse(str(input_path), "fasta"))
    queries: list[tuple[str, str]] = []
    for record in records:
        seq = str(record.seq).upper().replace("\n", "")
        seq = normalize_circular_query(reference_sequence, seq)
        queries.append((record.id, seq))

    aligned_pairs = align_batch_with_halign4(reference_sequence, queries, threads=threads)

    samples: list[SampleProfile] = []
    for record in records:
        aligned_reference, aligned_query = aligned_pairs[record.id]
        variants, ignored_positions, range_text = variants_from_alignment(aligned_reference, aligned_query)
        variants = apply_alignment_rules(variants, rules)
        samples.append(
            SampleProfile(
                sample_id=record.id,
                range_text=range_text,
                variants=ordered_unique([clean_token(token) for token in variants]),
                source_path=input_path,
                ignored_positions=ignored_positions,
            )
        )
    return samples
