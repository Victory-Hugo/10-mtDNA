"""FASTA 对齐与变异提取。"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner


_PROJECT_DIR = Path(__file__).resolve().parents[2]
HALIGN4_BIN = _PROJECT_DIR / "bin" / "halign4"
_HALIGN4_REF_ID = "__REFERENCE__"


def normalize_circular_query(reference: str, query: str) -> str:
    """把环形 mtDNA 序列旋转到与参考序列同起点。"""

    if not query or abs(len(reference) - len(query)) > 200:
        return query
    doubled_query = query + query
    for seed_length in (40, 30, 25, 20, 15, 12, 10):
        seed = reference[:seed_length]
        index = doubled_query.find(seed)
        if 0 <= index < len(query):
            return query[index:] + query[:index]
    return query


def global_align_to_reference(reference: str, query: str) -> tuple[str, str]:
    """使用 PairwiseAligner 做全局比对。"""

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -8.0
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(reference, query)
    if len(alignments) == 0:
        raise ValueError("FASTA 序列无法与参考序列完成全局比对。")
    return build_gapped_sequences(reference, query, alignments[0])


def build_gapped_sequences(reference: str, query: str, alignment) -> tuple[str, str]:
    """把 PairwiseAligner 的坐标片段重建成带 gap 的双序列。"""

    ref_parts: list[str] = []
    query_parts: list[str] = []
    ref_cursor = 0
    query_cursor = 0

    for (ref_start, ref_end), (query_start, query_end) in zip(alignment.aligned[0], alignment.aligned[1]):
        if ref_cursor < ref_start:
            ref_parts.append(reference[ref_cursor:ref_start])
            query_parts.append("-" * (ref_start - ref_cursor))
        if query_cursor < query_start:
            ref_parts.append("-" * (query_start - query_cursor))
            query_parts.append(query[query_cursor:query_start])
        ref_parts.append(reference[ref_start:ref_end])
        query_parts.append(query[query_start:query_end])
        ref_cursor = ref_end
        query_cursor = query_end

    if ref_cursor < len(reference):
        ref_parts.append(reference[ref_cursor:])
        query_parts.append("-" * (len(reference) - ref_cursor))
    if query_cursor < len(query):
        ref_parts.append("-" * (len(query) - query_cursor))
        query_parts.append(query[query_cursor:])

    return "".join(ref_parts), "".join(query_parts)


def align_batch_with_halign4(
    reference: str,
    queries: list[tuple[str, str]],
    threads: int = 1,
) -> dict[str, tuple[str, str]]:
    """用 halign4 批量比对所有查询序列到参考序列，返回 {sample_id: (aligned_ref, aligned_query)}。"""

    id_map: dict[str, str] = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = Path(tmpdir) / "input.fasta"
        output_fasta = Path(tmpdir) / "output.fasta"

        with input_fasta.open("w") as fh:
            fh.write(f">{_HALIGN4_REF_ID}\n{reference}\n")
            for idx, (sample_id, seq) in enumerate(queries):
                safe_id = f"seq_{idx}"
                id_map[safe_id] = sample_id
                fh.write(f">{safe_id}\n{seq}\n")

        subprocess.run(
            [str(HALIGN4_BIN), str(input_fasta), str(output_fasta),
             "-r", _HALIGN4_REF_ID, "-t", str(threads)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        aligned: dict[str, str] = {}
        for record in SeqIO.parse(str(output_fasta), "fasta"):
            aligned[record.id] = str(record.seq).upper()

    aligned_ref = aligned.get(_HALIGN4_REF_ID)
    if aligned_ref is None:
        raise ValueError("halign4 输出中未找到参考序列")

    return {
        id_map[safe_id]: (aligned_ref, aligned_query)
        for safe_id, aligned_query in aligned.items()
        if safe_id != _HALIGN4_REF_ID
    }


def variants_from_alignment(reference_aligned: str, query_aligned: str) -> tuple[list[str], set[int], str]:
    """从比对结果生成 HaploGrep 风格变异列表。"""

    variants: list[str] = []
    ignored_positions: set[int] = set()
    ref_position = 0
    insertion_buffer: list[str] = []
    insertion_anchor = 0
    covered_positions: list[int] = []

    for ref_base, query_base in zip(reference_aligned, query_aligned):
        if ref_base != "-":
            if insertion_buffer:
                variants.append(f"{insertion_anchor}.1{''.join(insertion_buffer)}")
                insertion_buffer = []
            ref_position += 1
            if query_base == "-":
                variants.append(f"{ref_position}d")
                continue
            if query_base.upper() == "N":
                ignored_positions.add(ref_position)
                continue
            covered_positions.append(ref_position)
            if ref_base.upper() != query_base.upper():
                variants.append(f"{ref_position}{query_base.upper()}")
            insertion_anchor = ref_position
            continue

        if query_base == "-":
            continue
        insertion_buffer.append(query_base.upper())

    if insertion_buffer:
        variants.append(f"{insertion_anchor}.1{''.join(insertion_buffer)}")

    if covered_positions:
        range_text = f"{covered_positions[0]}-{covered_positions[-1]}"
    else:
        range_text = "1-16569"
    return variants, ignored_positions, range_text
