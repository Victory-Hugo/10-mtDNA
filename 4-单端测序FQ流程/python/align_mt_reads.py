#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path
from typing import Sequence

log = logging.getLogger(__name__)


def run_command(cmd: Sequence[str]) -> None:
    log.info("Running command: %s", " ".join(cmd))
    subprocess.run(list(cmd), check=True)


def run(
    sample_id: str,
    fastq: str,
    reference: str,
    output_dir: str,
    threads: int,
    bwa_bin: str,
    samtools_bin: str,
    min_mapq: int,
) -> dict[str, str]:
    fastq_path = Path(fastq)
    reference_path = Path(reference)
    output_path = Path(output_dir)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ not found: {fastq_path}")
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference not found: {reference_path}")

    output_path.mkdir(parents=True, exist_ok=True)
    sorted_bam = output_path / f"{sample_id}.sorted.bam"
    flagstat_txt = output_path / f"{sample_id}.flagstat.txt"
    idxstats_txt = output_path / f"{sample_id}.idxstats.txt"
    primary_bam = output_path / f"{sample_id}.primary.q{min_mapq}.bam"
    primary_flagstat_txt = output_path / f"{sample_id}.primary.q{min_mapq}.flagstat.txt"
    primary_coverage_tsv = output_path / f"{sample_id}.primary.q{min_mapq}.coverage.tsv"

    bwa_cmd = [
        bwa_bin,
        "mem",
        "-t",
        str(threads),
        "-M",
        "-R",
        f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1",
        str(reference_path),
        str(fastq_path),
    ]
    sort_cmd = [
        samtools_bin,
        "sort",
        "-@",
        str(threads),
        "-o",
        str(sorted_bam),
        "-",
    ]
    log.info("Aligning sample %s", sample_id)
    bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
    try:
        sort_proc = subprocess.Popen(sort_cmd, stdin=bwa_proc.stdout)
        assert bwa_proc.stdout is not None
        bwa_proc.stdout.close()
        sort_rc = sort_proc.wait()
        bwa_rc = bwa_proc.wait()
    finally:
        if bwa_proc.stdout is not None:
            bwa_proc.stdout.close()
    if bwa_rc != 0:
        raise subprocess.CalledProcessError(bwa_rc, bwa_cmd)
    if sort_rc != 0:
        raise subprocess.CalledProcessError(sort_rc, sort_cmd)

    run_command([samtools_bin, "index", str(sorted_bam)])
    with open(flagstat_txt, "w", encoding="utf-8") as handle:
        subprocess.run([samtools_bin, "flagstat", str(sorted_bam)], check=True, stdout=handle)
    with open(idxstats_txt, "w", encoding="utf-8") as handle:
        subprocess.run([samtools_bin, "idxstats", str(sorted_bam)], check=True, stdout=handle)

    run_command(
        [
            samtools_bin,
            "view",
            "-@",
            str(threads),
            "-b",
            "-F",
            "2308",
            "-q",
            str(min_mapq),
            str(sorted_bam),
            "-o",
            str(primary_bam),
        ]
    )
    run_command([samtools_bin, "index", str(primary_bam)])
    with open(primary_flagstat_txt, "w", encoding="utf-8") as handle:
        subprocess.run([samtools_bin, "flagstat", str(primary_bam)], check=True, stdout=handle)
    with open(primary_coverage_tsv, "w", encoding="utf-8") as handle:
        subprocess.run([samtools_bin, "coverage", str(primary_bam)], check=True, stdout=handle)

    return {
        "sorted_bam": str(sorted_bam),
        "primary_bam": str(primary_bam),
        "flagstat": str(flagstat_txt),
        "primary_flagstat": str(primary_flagstat_txt),
        "primary_coverage": str(primary_coverage_tsv),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Align single-end mtDNA FASTQ reads to rCRS.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--bwa-bin", required=True)
    parser.add_argument("--samtools-bin", required=True)
    parser.add_argument("--min-mapq", type=int, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
