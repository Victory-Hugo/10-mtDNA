"""结果输出。"""

from __future__ import annotations

from pathlib import Path

from .models import ClassificationHit, SampleProfile, TreeBundle


def write_classification_report(
    output_path: Path,
    all_hits: dict[str, list[ClassificationHit]],
) -> None:
    """写默认 TSV 输出。"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write('"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\n')
        for sample_id, hits in all_hits.items():
            for hit in hits:
                handle.write(
                    f'"{sample_id}"\t"{hit.haplogroup}"\t"{hit.rank}"\t"{hit.quality:.4f}"\t"{hit.range_text}"\n'
                )


def write_extended_report(
    output_path: Path,
    all_hits: dict[str, list[ClassificationHit]],
) -> None:
    """写扩展 TSV 输出。"""

    extended_path = output_path.with_name(f"{output_path.stem}.extended{output_path.suffix or '.txt'}")
    with extended_path.open("w", encoding="utf-8") as handle:
        header = [
            "SampleID",
            "Haplogroup",
            "Rank",
            "Quality",
            "Range",
            "Input_Sample",
            "Found_Polymorphisms",
            "Missing_Polymorphisms",
            "Extra_Polymorphisms",
            "Found_Weight",
            "Expected_Weight",
            "Sample_Weight",
        ]
        handle.write("\t".join(f'"{column}"' for column in header) + "\n")
        for sample_id, hits in all_hits.items():
            for hit in hits:
                values = [
                    sample_id,
                    hit.haplogroup,
                    str(hit.rank),
                    f"{hit.quality:.4f}",
                    hit.range_text,
                    " ".join(sorted(set(hit.found_variants + hit.extra_variants))),
                    " ".join(hit.found_variants),
                    " ".join(hit.missing_variants),
                    " ".join(hit.extra_variants),
                    f"{hit.found_weight:.4f}",
                    f"{hit.expected_weight:.4f}",
                    f"{hit.sample_weight:.4f}",
                ]
                handle.write("\t".join(f'"{value}"' for value in values) + "\n")


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
