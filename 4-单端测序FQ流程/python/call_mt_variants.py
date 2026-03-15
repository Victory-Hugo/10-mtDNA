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
    input_bam: str,
    reference: str,
    rename_map: str,
    output_dir: str,
    threads: int,
    bcftools_bin: str,
    min_mapq: int,
    min_baseq: int,
) -> dict[str, str]:
    bam_path = Path(input_bam)
    reference_path = Path(reference)
    rename_map_path = Path(rename_map)
    output_path = Path(output_dir)
    for path in (bam_path, reference_path, rename_map_path):
        if not path.exists():
            raise FileNotFoundError(f"Required file not found: {path}")

    output_path.mkdir(parents=True, exist_ok=True)
    raw_chr_vcfgz = output_path / f"{sample_id}.raw.chrM.vcf.gz"
    raw_mt_vcf = output_path / f"{sample_id}.raw.MT.vcf"
    raw_stats_txt = output_path / f"{sample_id}.raw.stats.txt"

    mpileup_cmd = [
        bcftools_bin,
        "mpileup",
        "-f",
        str(reference_path),
        "-a",
        "FORMAT/DP,FORMAT/AD",
        "-q",
        str(min_mapq),
        "-Q",
        str(min_baseq),
        "-Ou",
        str(bam_path),
    ]
    call_cmd = [bcftools_bin, "call", "-mv", "--ploidy", "1", "-Ou"]
    norm_cmd = [
        bcftools_bin,
        "norm",
        "-f",
        str(reference_path),
        "-m",
        "-both",
        "-Oz",
        "-o",
        str(raw_chr_vcfgz),
    ]
    log.info("Calling variants for sample %s", sample_id)
    mpileup_proc = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
    try:
        call_proc = subprocess.Popen(call_cmd, stdin=mpileup_proc.stdout, stdout=subprocess.PIPE)
        assert mpileup_proc.stdout is not None
        mpileup_proc.stdout.close()
        norm_proc = subprocess.Popen(norm_cmd, stdin=call_proc.stdout)
        assert call_proc.stdout is not None
        call_proc.stdout.close()
        norm_rc = norm_proc.wait()
        call_rc = call_proc.wait()
        mpileup_rc = mpileup_proc.wait()
    finally:
        if mpileup_proc.stdout is not None:
            mpileup_proc.stdout.close()
        if call_proc.stdout is not None:
            call_proc.stdout.close()
    if mpileup_rc != 0:
        raise subprocess.CalledProcessError(mpileup_rc, mpileup_cmd)
    if call_rc != 0:
        raise subprocess.CalledProcessError(call_rc, call_cmd)
    if norm_rc != 0:
        raise subprocess.CalledProcessError(norm_rc, norm_cmd)

    run_command([bcftools_bin, "index", "-f", str(raw_chr_vcfgz)])
    run_command(
        [
            bcftools_bin,
            "annotate",
            "--rename-chrs",
            str(rename_map_path),
            str(raw_chr_vcfgz),
            "-Ov",
            "-o",
            str(raw_mt_vcf),
        ]
    )
    with open(raw_stats_txt, "w", encoding="utf-8") as handle:
        subprocess.run([bcftools_bin, "stats", str(raw_chr_vcfgz)], check=True, stdout=handle)

    return {
        "raw_chr_vcfgz": str(raw_chr_vcfgz),
        "raw_mt_vcf": str(raw_mt_vcf),
        "raw_stats": str(raw_stats_txt),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Call mtDNA variants from a filtered BAM.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--rename-map", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--bcftools-bin", required=True)
    parser.add_argument("--min-mapq", type=int, required=True)
    parser.add_argument("--min-baseq", type=int, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    run(**vars(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
