#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import statistics
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable


ADAPTER_MOTIFS = (
    "AGATCGGAAGAGC",
    "GCTCTTCCGATCT",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quick QC diagnosis for mtDNA FASTQ input."
    )
    parser.add_argument("--fastq", required=True, help="Input FASTQ(.gz) file.")
    parser.add_argument(
        "--outdir",
        default="/mnt/l/fastq/output",
        help="Root output directory.",
    )
    parser.add_argument(
        "--reference",
        default="/mnt/l/fastq/data/reference/chrM_rCRS.shifted.fa",
        help="Shifted chrM reference used for quick alignment checks.",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--dup-sample-size",
        type=int,
        default=200000,
        help="Number of reads used for exact-duplicate estimation.",
    )
    parser.add_argument(
        "--sample-name",
        default=None,
        help="Optional sample name override.",
    )
    return parser.parse_args()


def derive_sample_name(fastq: Path, override: str | None) -> str:
    if override:
        return override
    name = fastq.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return fastq.stem


def run_command(
    cmd: list[str],
    *,
    stdout_path: Path | None = None,
    stderr_path: Path | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    stdout_handle = open(stdout_path, "w", encoding="utf-8") if stdout_path else subprocess.PIPE
    stderr_handle = open(stderr_path, "w", encoding="utf-8") if stderr_path else subprocess.PIPE
    try:
        result = subprocess.run(
            cmd,
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=check,
        )
    finally:
        if stdout_path:
            stdout_handle.close()
        if stderr_path:
            stderr_handle.close()
    return result


def percentile(sorted_values: list[float], pct: float) -> float:
    if not sorted_values:
        return 0.0
    rank = math.ceil(len(sorted_values) * pct / 100.0) - 1
    rank = max(0, min(rank, len(sorted_values) - 1))
    return sorted_values[rank]


def open_fastq(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def analyze_fastq(fastq: Path, dup_sample_size: int) -> tuple[dict[str, object], list[dict[str, object]]]:
    read_count = 0
    total_bases = 0
    lengths: list[int] = []
    mean_qualities: list[float] = []
    base_counts: Counter[str] = Counter()
    adapter_counts: Counter[str] = Counter()
    exact_seq_counts: Counter[str] = Counter()
    poly_a_reads = 0
    poly_t_reads = 0
    reads_lt_100 = 0
    reads_lt_140 = 0
    q20_bases = 0
    q30_bases = 0
    n_bases = 0
    top_headers: list[str] = []

    with open_fastq(fastq) as handle:
        while True:
            header = handle.readline().rstrip()
            if not header:
                break
            seq = handle.readline().rstrip()
            handle.readline()
            qual = handle.readline().rstrip()

            read_count += 1
            if len(top_headers) < 5:
                top_headers.append(header)

            read_len = len(seq)
            lengths.append(read_len)
            total_bases += read_len
            base_counts.update(seq)
            n_bases += seq.count("N")

            if read_len < 100:
                reads_lt_100 += 1
            if read_len < 140:
                reads_lt_140 += 1

            quals = [ord(ch) - 33 for ch in qual]
            mean_q = sum(quals) / read_len if read_len else 0.0
            mean_qualities.append(mean_q)
            q20_bases += sum(1 for q in quals if q >= 20)
            q30_bases += sum(1 for q in quals if q >= 30)

            for motif in ADAPTER_MOTIFS:
                if motif in seq:
                    adapter_counts[motif] += 1

            if "AAAAAAAAAA" in seq:
                poly_a_reads += 1
            if "TTTTTTTTTT" in seq:
                poly_t_reads += 1

            if read_count <= dup_sample_size:
                exact_seq_counts[seq] += 1

    if read_count == 0:
        raise ValueError(f"No reads found in {fastq}")

    length_counts = Counter(lengths)
    sorted_mean_q = sorted(mean_qualities)
    base_total = sum(base_counts.values())
    dup_total = sum(exact_seq_counts.values())
    dup_unique = len(exact_seq_counts)
    top_sequences = []

    for seq, count in exact_seq_counts.most_common(20):
        top_sequences.append(
            {
                "sequence": seq,
                "count": count,
                "fraction_in_duplicate_sample": round(count / dup_total, 6) if dup_total else 0.0,
                "contains_adapter": any(motif in seq for motif in ADAPTER_MOTIFS),
            }
        )

    metrics: dict[str, object] = {
        "read_count": read_count,
        "total_bases": total_bases,
        "min_length": min(lengths),
        "mean_length": round(statistics.mean(lengths), 3),
        "median_length": round(statistics.median(lengths), 3),
        "max_length": max(lengths),
        "mode_length": length_counts.most_common(1)[0][0],
        "reads_lt_100": reads_lt_100,
        "reads_lt_140": reads_lt_140,
        "reads_lt_100_fraction": round(reads_lt_100 / read_count, 6),
        "reads_lt_140_fraction": round(reads_lt_140 / read_count, 6),
        "mean_read_quality": round(statistics.mean(mean_qualities), 3),
        "median_read_quality": round(statistics.median(mean_qualities), 3),
        "read_quality_p01": round(percentile(sorted_mean_q, 1), 3),
        "read_quality_p05": round(percentile(sorted_mean_q, 5), 3),
        "read_quality_p10": round(percentile(sorted_mean_q, 10), 3),
        "read_quality_p25": round(percentile(sorted_mean_q, 25), 3),
        "read_quality_p50": round(percentile(sorted_mean_q, 50), 3),
        "read_quality_p75": round(percentile(sorted_mean_q, 75), 3),
        "read_quality_p90": round(percentile(sorted_mean_q, 90), 3),
        "read_quality_p95": round(percentile(sorted_mean_q, 95), 3),
        "read_quality_p99": round(percentile(sorted_mean_q, 99), 3),
        "q20_base_fraction": round(q20_bases / total_bases, 6),
        "q30_base_fraction": round(q30_bases / total_bases, 6),
        "base_fraction_A": round(base_counts["A"] / base_total, 6),
        "base_fraction_C": round(base_counts["C"] / base_total, 6),
        "base_fraction_G": round(base_counts["G"] / base_total, 6),
        "base_fraction_T": round(base_counts["T"] / base_total, 6),
        "base_fraction_N": round(base_counts["N"] / base_total, 6),
        "adapter_reads_forward": adapter_counts[ADAPTER_MOTIFS[0]],
        "adapter_reads_reverse": adapter_counts[ADAPTER_MOTIFS[1]],
        "adapter_reads_forward_fraction": round(adapter_counts[ADAPTER_MOTIFS[0]] / read_count, 6),
        "adapter_reads_reverse_fraction": round(adapter_counts[ADAPTER_MOTIFS[1]] / read_count, 6),
        "polyA10_reads": poly_a_reads,
        "polyA10_fraction": round(poly_a_reads / read_count, 6),
        "polyT10_reads": poly_t_reads,
        "polyT10_fraction": round(poly_t_reads / read_count, 6),
        "dup_sample_reads": dup_total,
        "dup_sample_unique_reads": dup_unique,
        "dup_sample_exact_duplicate_fraction": round((dup_total - dup_unique) / dup_total, 6) if dup_total else 0.0,
        "first_headers": top_headers,
    }
    return metrics, top_sequences


def write_key_value_table(path: Path, metrics: dict[str, object]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, json.dumps(value, ensure_ascii=False) if isinstance(value, list) else value])


def write_top_sequences(path: Path, rows: list[dict[str, object]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sequence", "count", "fraction_in_duplicate_sample", "contains_adapter"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def run_seqkit_outputs(fastq: Path, outdir: Path) -> None:
    run_command(
        ["seqkit", "stats", "-Ta", str(fastq)],
        stdout_path=outdir / "seqkit.stats.tsv",
    )
    fqchk_path = outdir / "seqkit.fqchk.txt"
    try:
        run_command(
            ["seqkit", "fqchk", "-q", "33", str(fastq)],
            stdout_path=fqchk_path,
        )
    except subprocess.CalledProcessError as exc:
        with open(fqchk_path, "w", encoding="utf-8") as handle:
            handle.write(
                "seqkit fqchk is unavailable in the current seqkit build; "
                "basic FASTQ metrics were still collected by the Python analyzer.\n"
            )
            if exc.stderr:
                handle.write(str(exc.stderr))


def run_alignment(fastq: Path, reference: Path, outdir: Path, threads: int) -> Path:
    bam_path = outdir / "chrM.quickcheck.sorted.bam"
    bwa_log = outdir / "chrM.quickcheck.bwa.log"
    with open(bam_path, "wb") as bam_handle, open(bwa_log, "w", encoding="utf-8") as log_handle:
        bwa = subprocess.Popen(
            ["bwa", "mem", "-t", str(threads), str(reference), str(fastq)],
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        sort_proc = subprocess.Popen(
            ["samtools", "sort", "-@", str(threads), "-o", str(bam_path), "-"],
            stdin=bwa.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if bwa.stdout is not None:
            bwa.stdout.close()
        _, sort_stderr = sort_proc.communicate()
        bwa_rc = bwa.wait()
        if bwa_rc != 0:
            raise subprocess.CalledProcessError(bwa_rc, bwa.args)
        if sort_proc.returncode != 0:
            raise subprocess.CalledProcessError(
                sort_proc.returncode,
                sort_proc.args,
                stderr=sort_stderr.decode("utf-8", errors="replace"),
            )
    run_command(["samtools", "index", str(bam_path)])
    run_command(["samtools", "flagstat", str(bam_path)], stdout_path=outdir / "chrM.quickcheck.flagstat.txt")
    run_command(["samtools", "coverage", str(bam_path)], stdout_path=outdir / "chrM.quickcheck.coverage.tsv")
    run_command(["samtools", "idxstats", str(bam_path)], stdout_path=outdir / "chrM.quickcheck.idxstats.tsv")
    return bam_path


def parse_flagstat(path: Path) -> dict[str, int]:
    metrics: dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            match = re.match(r"^(\d+) \+ \d+ (.+)$", line.strip())
            if not match:
                continue
            value = int(match.group(1))
            label = match.group(2)
            if label.startswith("in total"):
                metrics["alignments_total"] = value
            elif label.startswith("primary mapped"):
                metrics["primary_mapped_alignments"] = value
            elif label == "primary":
                metrics["primary_alignments"] = value
            elif label.startswith("secondary"):
                metrics["secondary_alignments"] = value
            elif label.startswith("supplementary"):
                metrics["supplementary_alignments"] = value
            elif label.startswith("mapped"):
                metrics["mapped_alignments"] = value
            elif label.startswith("paired in sequencing"):
                metrics["paired_in_sequencing"] = value
    return metrics


def parse_coverage(path: Path) -> dict[str, float]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        return {}
    row = rows[0]
    return {
        "coverage_percent": float(row["coverage"]),
        "mean_depth": float(row["meandepth"]),
        "mean_baseq": float(row["meanbaseq"]),
        "mean_mapq": float(row["meanmapq"]),
        "covered_bases": float(row["covbases"]),
        "reference_start": float(row["startpos"]),
        "reference_end": float(row["endpos"]),
        "reads_on_reference": float(row["numreads"]),
    }


def summarize_primary_alignments(bam_path: Path, outdir: Path) -> dict[str, object]:
    primary_reads = 0
    primary_mapped = 0
    primary_clip_ge20 = 0
    primary_with_sa = 0
    supplementary_examples: list[dict[str, object]] = []
    primary_examples: list[dict[str, object]] = []

    proc = subprocess.Popen(
        ["samtools", "view", str(bam_path)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        qname = fields[0]
        flag = int(fields[1])
        rname = fields[2]
        pos = int(fields[3])
        cigar = fields[5]
        sa_tag = next((field for field in fields[11:] if field.startswith("SA:Z:")), "")

        if flag & 2048:
            if len(supplementary_examples) < 10:
                supplementary_examples.append(
                    {
                        "category": "supplementary",
                        "qname": qname,
                        "flag": flag,
                        "rname": rname,
                        "pos": pos,
                        "cigar": cigar,
                        "sa_tag": sa_tag,
                    }
                )
            continue

        primary_reads += 1
        if cigar != "*":
            primary_mapped += 1

        clip_lengths = [int(match[:-1]) for match in re.findall(r"\d+[SH]", cigar)]
        has_big_clip = any(length >= 20 for length in clip_lengths)
        if has_big_clip:
            primary_clip_ge20 += 1
            if len(primary_examples) < 10:
                primary_examples.append(
                    {
                        "category": "primary_clipped",
                        "qname": qname,
                        "flag": flag,
                        "rname": rname,
                        "pos": pos,
                        "cigar": cigar,
                        "sa_tag": sa_tag,
                    }
                )
        if sa_tag:
            primary_with_sa += 1

    rc = proc.wait()
    if rc != 0:
        raise subprocess.CalledProcessError(rc, proc.args)

    summary = {
        "primary_reads": primary_reads,
        "primary_mapped": primary_mapped,
        "primary_clip_ge20": primary_clip_ge20,
        "primary_clip_ge20_fraction": round(primary_clip_ge20 / primary_reads, 6) if primary_reads else 0.0,
        "primary_with_sa": primary_with_sa,
        "primary_with_sa_fraction": round(primary_with_sa / primary_reads, 6) if primary_reads else 0.0,
    }

    write_key_value_table(outdir / "chrM.quickcheck.primary_summary.tsv", summary)
    with open(outdir / "chrM.quickcheck.alignment_examples.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["category", "qname", "flag", "rname", "pos", "cigar", "sa_tag"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(primary_examples + supplementary_examples)
    return summary


def build_report(
    sample_name: str,
    fastq_path: Path,
    outdir: Path,
    fastq_metrics: dict[str, object],
    top_sequences: list[dict[str, object]],
    flagstat_metrics: dict[str, int],
    coverage_metrics: dict[str, float],
    primary_summary: dict[str, object],
) -> str:
    adapter_fraction = float(fastq_metrics["adapter_reads_forward_fraction"])
    q30_fraction = float(fastq_metrics["q30_base_fraction"])
    primary_map_fraction = (
        flagstat_metrics.get("primary_mapped_alignments", 0) / int(fastq_metrics["read_count"])
        if fastq_metrics["read_count"]
        else 0.0
    )
    supplementary_fraction = (
        flagstat_metrics.get("supplementary_alignments", 0) / flagstat_metrics.get("alignments_total", 1)
        if flagstat_metrics.get("alignments_total")
        else 0.0
    )
    clipped_fraction = float(primary_summary["primary_clip_ge20_fraction"])
    sa_fraction = float(primary_summary["primary_with_sa_fraction"])

    top_adapter_seq = next((row for row in top_sequences if row["contains_adapter"]), None)
    adapter_line = (
        f"The duplicate sample includes an overrepresented adapter-like sequence seen {top_adapter_seq['count']} times."
        if top_adapter_seq
        else "No adapter-like sequence dominated the duplicate sample."
    )

    issues: list[str] = []
    if adapter_fraction >= 0.01:
        issues.append(
            f"Adapter contamination is non-trivial: {adapter_fraction:.2%} of reads contain the forward Illumina adapter motif."
        )
    if q30_fraction < 0.6:
        issues.append(
            f"Per-base quality is moderate rather than excellent: Q30 bases are {q30_fraction:.2%}."
        )
    if primary_map_fraction < 0.5:
        issues.append(
            f"Only {primary_map_fraction:.2%} of reads produce a primary chrM alignment."
        )
    if clipped_fraction >= 0.2 or sa_fraction >= 0.2:
        issues.append(
            "A large fraction of primary alignments are heavily clipped or have supplementary alignments, consistent with circular-boundary spanning reads or concatemeric/chimeric fragments."
        )
    if not issues:
        issues.append("No obvious failure mode was detected in the raw read statistics.")

    lines = [
        f"# Quick mtDNA FASTQ diagnosis: {sample_name}",
        "",
        "## Input",
        f"- FASTQ: `{fastq_path}`",
        f"- Output directory: `{outdir}`",
        "- Execution mode: single FASTQ quick diagnosis, not a full pipeline replacement.",
        "",
        "## Read-level summary",
        f"- Reads: {int(fastq_metrics['read_count']):,}",
        f"- Length: min {int(fastq_metrics['min_length'])}, mean {float(fastq_metrics['mean_length']):.2f}, median {float(fastq_metrics['median_length']):.2f}, max {int(fastq_metrics['max_length'])}",
        f"- Mean read quality: {float(fastq_metrics['mean_read_quality']):.2f}",
        f"- Q20 bases: {float(fastq_metrics['q20_base_fraction']):.2%}",
        f"- Q30 bases: {float(fastq_metrics['q30_base_fraction']):.2%}",
        f"- Reads <100 bp: {int(fastq_metrics['reads_lt_100']):,} ({float(fastq_metrics['reads_lt_100_fraction']):.2%})",
        "",
        "## Contamination / structure signals",
        f"- Reads containing `AGATCGGAAGAGC`: {int(fastq_metrics['adapter_reads_forward']):,} ({float(fastq_metrics['adapter_reads_forward_fraction']):.2%})",
        f"- Reads containing `GCTCTTCCGATCT`: {int(fastq_metrics['adapter_reads_reverse']):,} ({float(fastq_metrics['adapter_reads_reverse_fraction']):.2%})",
        f"- Reads with polyA(10): {int(fastq_metrics['polyA10_reads']):,} ({float(fastq_metrics['polyA10_fraction']):.2%})",
        f"- Reads with polyT(10): {int(fastq_metrics['polyT10_reads']):,} ({float(fastq_metrics['polyT10_fraction']):.2%})",
        f"- Exact duplicate fraction in the first {int(fastq_metrics['dup_sample_reads']):,} reads: {float(fastq_metrics['dup_sample_exact_duplicate_fraction']):.2%}",
        f"- {adapter_line}",
        "",
        "## Quick chrM alignment check",
        f"- Total alignments: {flagstat_metrics.get('alignments_total', 0):,}",
        f"- Primary alignments: {flagstat_metrics.get('primary_alignments', 0):,}",
        f"- Primary mapped alignments: {flagstat_metrics.get('primary_mapped_alignments', 0):,} ({primary_map_fraction:.2%} of reads)",
        f"- Supplementary alignments: {flagstat_metrics.get('supplementary_alignments', 0):,} ({supplementary_fraction:.2%} of all alignments)",
        f"- `paired in sequencing` from BAM flags: {flagstat_metrics.get('paired_in_sequencing', 0):,}",
        f"- chrM coverage breadth: {coverage_metrics.get('coverage_percent', 0.0):.2f}%",
        f"- chrM mean depth: {coverage_metrics.get('mean_depth', 0.0):.2f}x",
        f"- Primary reads with >=20 bp clipping: {int(primary_summary['primary_clip_ge20']):,} ({clipped_fraction:.2%})",
        f"- Primary reads carrying an `SA` tag: {int(primary_summary['primary_with_sa']):,} ({sa_fraction:.2%})",
        "",
        "## Interpretation",
    ]
    lines.extend(f"- {issue}" for issue in issues)
    lines.extend(
        [
            "- chrM signal is present and deep, so this sample is not simply an empty or failed mtDNA library.",
            "- The current `4-Mutect2` workflow expects BAM input and later reconstructs `R1/R2` FASTQs for paired-end realignment. A single FASTQ input like this should be treated as a workflow-compatibility risk unless an actual mate file exists upstream.",
            "- The clipped/supplementary alignment pattern around chrM boundaries can inflate artifacts in variant calling if handled as standard genomic short reads.",
            "",
            "## Generated files",
            "- `fastq_metrics.tsv`: read-level metrics computed directly from the FASTQ",
            "- `top_overrepresented_sequences.tsv`: the most frequent exact sequences in the duplicate sample",
            "- `seqkit.stats.tsv` and `seqkit.fqchk.txt`: raw seqkit QC outputs",
            "- `chrM.quickcheck.*`: alignment BAM, index, flagstat, coverage, idxstats, and example alignments",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    fastq = Path(args.fastq).resolve()
    reference = Path(args.reference).resolve()
    if not fastq.exists():
        raise FileNotFoundError(f"FASTQ not found: {fastq}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference not found: {reference}")

    sample_name = derive_sample_name(fastq, args.sample_name)
    sample_outdir = Path(args.outdir).resolve() / sample_name
    sample_outdir.mkdir(parents=True, exist_ok=True)

    fastq_metrics, top_sequences = analyze_fastq(fastq, args.dup_sample_size)
    write_key_value_table(sample_outdir / "fastq_metrics.tsv", fastq_metrics)
    write_top_sequences(sample_outdir / "top_overrepresented_sequences.tsv", top_sequences)
    with open(sample_outdir / "fastq_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(fastq_metrics, handle, ensure_ascii=False, indent=2)

    run_seqkit_outputs(fastq, sample_outdir)
    bam_path = run_alignment(fastq, reference, sample_outdir, args.threads)
    flagstat_metrics = parse_flagstat(sample_outdir / "chrM.quickcheck.flagstat.txt")
    coverage_metrics = parse_coverage(sample_outdir / "chrM.quickcheck.coverage.tsv")
    primary_summary = summarize_primary_alignments(bam_path, sample_outdir)

    report = build_report(
        sample_name,
        fastq,
        sample_outdir,
        fastq_metrics,
        top_sequences,
        flagstat_metrics,
        coverage_metrics,
        primary_summary,
    )
    with open(sample_outdir / "diagnosis_report.md", "w", encoding="utf-8") as handle:
        handle.write(report)
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
