"""结果输出。"""

from __future__ import annotations

from pathlib import Path

from .models import ClassificationHit, SampleProfile, TreeBundle

_BASIC_HEADER = '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\n'

_EXTENDED_COLUMNS = [
    "SampleID", "Haplogroup", "Rank", "Quality", "Range",
    "Input_Sample", "Found_Polymorphisms", "Missing_Polymorphisms",
    "Extra_Polymorphisms", "AAC_In_Remainings",
    "Found_Weight", "Expected_Weight", "Sample_Weight",
]
_EXTENDED_HEADER = "\t".join(f'"{c}"' for c in _EXTENDED_COLUMNS) + "\n"


def _extended_path(output_path: Path) -> Path:
    return output_path.with_name(f"{output_path.stem}.extended{output_path.suffix or '.txt'}")


def write_classification_report(
    output_path: Path,
    all_hits: dict[str, list[ClassificationHit]],
    append: bool = False,
) -> None:
    """写默认 TSV 输出。append=True 时跳过表头并以追加模式写入。"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("a" if append else "w", encoding="utf-8") as handle:
        if not append:
            handle.write(_BASIC_HEADER)
        for sample_id, hits in all_hits.items():
            for hit in hits:
                handle.write(
                    f'"{sample_id}"\t"{hit.haplogroup}"\t"{hit.rank}"\t"{hit.quality:.4f}"\t"{hit.range_text}"\n'
                )


def write_extended_report(
    output_path: Path,
    all_hits: dict[str, list[ClassificationHit]],
    append: bool = False,
) -> None:
    """写扩展 TSV 输出。append=True 时跳过表头并以追加模式写入。"""

    ext_path = _extended_path(output_path)
    with ext_path.open("a" if append else "w", encoding="utf-8") as handle:
        if not append:
            handle.write(_EXTENDED_HEADER)
        for sample_id, hits in all_hits.items():
            for hit in hits:
                all_found = sorted(
                    set(hit.found_variants) | set(hit.confirmed_back_mutations),
                    key=lambda t: (int("".join(c for c in t if c.isdigit()) or "0"), t),
                )
                annotated_extra = " ".join(
                    f"{v} ({ann})" for v, ann in hit.annotated_extra_variants
                )
                values = [
                    sample_id,
                    hit.haplogroup,
                    str(hit.rank),
                    f"{hit.quality:.4f}",
                    hit.range_text,
                    " ".join(sorted(set(hit.found_variants + hit.extra_variants))),
                    " ".join(all_found),
                    " ".join(hit.missing_variants),
                    annotated_extra,
                    " ".join(hit.aac_in_remainings),
                    f"{hit.found_weight:.4f}",
                    f"{hit.expected_weight:.4f}",
                    f"{hit.sample_weight:.4f}",
                ]
                handle.write("\t".join(f'"{v}"' for v in values) + "\n")


def write_reconstructed_fasta(
    output_path: Path,
    samples: list[SampleProfile],
    tree_bundle: TreeBundle,
) -> None:
    """根据样本变异重建序列。"""

    fasta_path = output_path.with_suffix(".fasta")
    with fasta_path.open("w", encoding="utf-8") as handle:
        for sample in samples:
            sequence = list(tree_bundle.reference_sequence)
            for token in sample.variants:
                if "." in token or token.endswith("d"):
                    continue
                digits = ""
                letters = ""
                for character in token:
                    if character.isdigit():
                        digits += character
                    else:
                        letters += character
                if not digits or not letters:
                    continue
                sequence[int(digits) - 1] = letters[-1]
            handle.write(f">{sample.sample_id}\n")
            for start in range(0, len(sequence), 70):
                handle.write("".join(sequence[start : start + 70]) + "\n")
