#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import logging
import os
import subprocess
import tempfile
import time
from typing import Dict, List, Optional, Tuple

import yaml


def load_config(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def setup_logger(log_path: str) -> logging.Logger:
    logger = logging.getLogger(log_path)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    if not logger.handlers:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        handler = logging.FileHandler(log_path)
        formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def log_status(logger: logging.Logger, status: str, reason: str = "") -> None:
    if reason:
        logger.info("STATUS\t%s\t%s", status, reason)
    else:
        logger.info("STATUS\t%s", status)


def run_command(command: List[str], logger: logging.Logger) -> None:
    logger.info("CMD\t%s", " ".join(command))
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.stdout:
        logger.info("STDOUT\t%s", result.stdout.strip())
    if result.stderr:
        logger.info("STDERR\t%s", result.stderr.strip())
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(command)}")


def get_samtools_version(samtools: str, logger: logging.Logger) -> None:
    try:
        result = subprocess.run([samtools, "--version"], check=False, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"samtools not found: {samtools}") from exc
    if result.stdout:
        first_line = result.stdout.splitlines()[0]
        logger.info("SAMTOOLS_VERSION\t%s", first_line)
    if result.returncode != 0:
        raise RuntimeError("samtools --version failed")


def get_mosdepth_version(mosdepth: str, logger: logging.Logger) -> None:
    try:
        result = subprocess.run([mosdepth, "--version"], check=False, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"mosdepth not found: {mosdepth}") from exc
    version_line = result.stdout.strip() or result.stderr.strip()
    if version_line:
        logger.info("MOSDEPTH_VERSION\t%s", version_line.splitlines()[0])
    if result.returncode != 0:
        raise RuntimeError("mosdepth --version failed")


def read_fai(fai_path: str) -> Dict[str, int]:
    contigs: Dict[str, int] = {}
    with open(fai_path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            contigs[parts[0]] = int(parts[1])
    return contigs


def detect_prefix(contigs: Dict[str, int]) -> str:
    if "chr1" in contigs:
        return "chr"
    return ""


def detect_mt_contig(contigs: Dict[str, int]) -> str:
    candidates = ["chrM", "MT", "chrMT", "M", "chrM"]
    for name in candidates:
        if name in contigs:
            return name
    raise RuntimeError("Cannot find mitochondrial contig in reference index")


def build_autosome_list(contigs: Dict[str, int], prefix: str) -> List[str]:
    autosomes = [f"{prefix}{i}" for i in range(1, 23)]
    missing = [name for name in autosomes if name not in contigs]
    if missing:
        raise RuntimeError(f"Missing autosomes in reference: {','.join(missing)}")
    return autosomes


def autosome_bed_path(tmp_dir: str, reference_fasta: str) -> str:
    key = hashlib.md5(reference_fasta.encode("utf-8")).hexdigest()[:10]
    return os.path.join(tmp_dir, f"autosomes_{key}.bed")


def ensure_autosome_bed(
    tmp_dir: str,
    contigs: Dict[str, int],
    autosomes: List[str],
    reference_fasta: str,
) -> str:
    os.makedirs(tmp_dir, exist_ok=True)
    bed_path = autosome_bed_path(tmp_dir, reference_fasta)
    if os.path.exists(bed_path):
        return bed_path

    for attempt in range(3):
        tmp_path = None
        try:
            os.makedirs(tmp_dir, exist_ok=True)
            fd, tmp_path = tempfile.mkstemp(prefix="autosomes_", suffix=".bed.tmp", dir=tmp_dir)
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                for name in autosomes:
                    handle.write(f"{name}\t0\t{contigs[name]}\n")
            os.replace(tmp_path, bed_path)
            return bed_path
        except FileNotFoundError:
            if attempt == 2:
                raise
            time.sleep(0.2)
        finally:
            if tmp_path and os.path.exists(tmp_path):
                os.remove(tmp_path)

    return bed_path


def run_mosdepth(command: List[str], logger: logging.Logger) -> None:
    logger.info("CMD\t%s", " ".join(command))
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.stdout:
        logger.info("STDOUT\t%s", result.stdout.strip())
    if result.stderr:
        logger.info("STDERR\t%s", result.stderr.strip())
    if result.returncode != 0:
        raise RuntimeError(f"mosdepth failed: {' '.join(command)}")


def read_mosdepth_summary(summary_path: str, contig: str) -> float:
    with open(summary_path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if not parts or parts[0] == "chrom":
                continue
            if parts[0] == contig:
                return float(parts[3])
    raise RuntimeError(f"Missing contig in mosdepth summary: {contig}")


def read_mosdepth_regions(regions_path: str) -> Tuple[float, int]:
    total_weighted = 0.0
    total_len = 0
    with gzip.open(regions_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            start = int(parts[1])
            end = int(parts[2])
            mean = float(parts[3])
            length = end - start
            total_weighted += mean * length
            total_len += length
    if total_len == 0:
        raise RuntimeError("No regions found in mosdepth output")
    return total_weighted / total_len, total_len


def mosdepth_prefix(tmp_dir: str, sample_id: str, tag: str) -> str:
    return os.path.join(tmp_dir, sample_id, f"{sample_id}.{tag}")


def should_skip(log_path: str) -> bool:
    if not os.path.exists(log_path):
        return False
    with open(log_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if "STATUS\tSUCCESS" in line:
                return True
    return False


def append_success_log(success_log: Optional[str], sample_id: str) -> None:
    if not success_log:
        return
    os.makedirs(os.path.dirname(success_log), exist_ok=True)
    with open(success_log, "a", encoding="utf-8") as handle:
        handle.write(f"{sample_id}\n")


def run(
    sample_id: str,
    bam_path: str,
    reference_fasta: str,
    output_path: str,
    log_path: str,
    tmp_dir: str,
    samtools: str = "samtools",
    mosdepth: str = "mosdepth",
    mosdepth_threads: int = 4,
    success_log: Optional[str] = None,
    force: bool = False,
) -> None:
    logger = setup_logger(log_path)
    if should_skip(log_path) and not force:
        log_status(logger, "SKIP", "previous success")
        return

    logger.info("SAMPLE\t%s", sample_id)
    logger.info("BAM\t%s", bam_path)
    logger.info("REFERENCE\t%s", reference_fasta)
    logger.info("OUTPUT\t%s", output_path)
    logger.info("TMP\t%s", tmp_dir)

    fai_path = f"{reference_fasta}.fai"
    if not os.path.exists(fai_path):
        run_command([samtools, "faidx", reference_fasta], logger)

    get_samtools_version(samtools, logger)
    get_mosdepth_version(mosdepth, logger)

    contigs = read_fai(fai_path)
    prefix = detect_prefix(contigs)
    autosomes = build_autosome_list(contigs, prefix)
    mt_contig = detect_mt_contig(contigs)
    bed_path = ensure_autosome_bed(tmp_dir, contigs, autosomes, reference_fasta)

    sample_tmp_dir = os.path.join(tmp_dir, sample_id)
    os.makedirs(sample_tmp_dir, exist_ok=True)

    mt_prefix = mosdepth_prefix(tmp_dir, sample_id, "mt")
    mt_command = [
        mosdepth,
        "--no-per-base",
        "--threads",
        str(mosdepth_threads),
        "--chrom",
        mt_contig,
        mt_prefix,
        bam_path,
    ]
    run_mosdepth(mt_command, logger)
    mt_summary = f"{mt_prefix}.mosdepth.summary.txt"
    mean_mt = read_mosdepth_summary(mt_summary, mt_contig)

    auto_prefix = mosdepth_prefix(tmp_dir, sample_id, "autosome")
    auto_command = [
        mosdepth,
        "--no-per-base",
        "--threads",
        str(mosdepth_threads),
        "--by",
        bed_path,
        auto_prefix,
        bam_path,
    ]
    run_mosdepth(auto_command, logger)
    regions_path = f"{auto_prefix}.regions.bed.gz"
    mean_auto, auto_count = read_mosdepth_regions(regions_path)
    logger.info("DEPTH_COUNT_AUTO\t%d", auto_count)

    mt_copy_number = (mean_mt / mean_auto) * 2

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_output = os.path.join(tmp_dir, f"{sample_id}.mt_copy_number.tmp")
    with open(tmp_output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "mean_mtDNA_depth",
                "mean_autosomal_depth",
                "mtDNA_copy_number",
            ]
        )
        writer.writerow(
            [
                sample_id,
                f"{mean_mt:.6f}",
                f"{mean_auto:.6f}",
                f"{mt_copy_number:.6f}",
            ]
        )
    os.replace(tmp_output, output_path)
    append_success_log(success_log, sample_id)
    log_status(logger, "SUCCESS")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute mtDNA copy number from a BAM file.")
    parser.add_argument("--config", help="Path to YAML config")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--bam", required=True, help="Input BAM path")
    parser.add_argument("--reference-fasta", help="Reference FASTA path")
    parser.add_argument("--output", help="Output TSV path")
    parser.add_argument("--log", help="Log file path")
    parser.add_argument("--tmp-dir", help="Temporary directory")
    parser.add_argument("--samtools", help="samtools executable")
    parser.add_argument("--mosdepth", help="mosdepth executable")
    parser.add_argument("--mosdepth-threads", type=int, help="mosdepth threads")
    parser.add_argument("--success-log", help="Path to success log")
    parser.add_argument("--force", action="store_true", help="Force recompute")
    args = parser.parse_args()

    config = load_config(args.config)
    reference_fasta = args.reference_fasta or config.get("reference_fasta")
    output_dir = config.get("output_dir")
    tmp_dir = args.tmp_dir or config.get("tmp_dir")
    log_dir = config.get("log_dir")
    samtools = args.samtools or config.get("samtools", "samtools")
    mosdepth = args.mosdepth or config.get("mosdepth", "mosdepth")
    mosdepth_threads = args.mosdepth_threads or int(config.get("mosdepth_threads", 4))
    success_log = args.success_log or config.get("success_log")

    if not reference_fasta:
        raise SystemExit("Missing reference_fasta")
    if not tmp_dir:
        raise SystemExit("Missing tmp_dir")
    if not output_dir:
        raise SystemExit("Missing output_dir")
    if not log_dir:
        raise SystemExit("Missing log_dir")

    output_path = args.output or os.path.join(output_dir, f"{args.sample_id}.tsv")
    log_path = args.log or os.path.join(log_dir, f"mt_copy_number_{args.sample_id}.log")

    try:
        run(
            sample_id=args.sample_id,
            bam_path=args.bam,
            reference_fasta=reference_fasta,
            output_path=output_path,
            log_path=log_path,
            tmp_dir=tmp_dir,
            samtools=samtools,
            mosdepth=mosdepth,
            mosdepth_threads=mosdepth_threads,
            success_log=success_log,
            force=args.force,
        )
    except Exception as exc:
        logger = setup_logger(log_path)
        log_status(logger, "FAILED", str(exc))
        raise


if __name__ == "__main__":
    main()
