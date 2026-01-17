#!/usr/bin/env python3
import argparse
import datetime as dt
import subprocess
from pathlib import Path

import yaml


FREEMIX_FAIL_THRESHOLD = 0.02
TARGET_CONTIGS = {str(i) for i in range(1, 23)} | {"X", "Y", "M"}


def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    if "verifybamid2" not in cfg:
        raise ValueError("Missing 'verifybamid2' section in config")
    return cfg["verifybamid2"]


def log_has_success(log_path: Path) -> bool:
    if not log_path.exists():
        return False
    try:
        content = log_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return False
    return "STATUS=SUCCESS" in content


def build_command(cfg: dict, bam: str, output_prefix: str) -> list:
    cmd = [
        "conda",
        "run",
        "-n",
        cfg["conda_env"],
        cfg["binary"],
        "--BamFile",
        bam,
        "--Reference",
        cfg["reference"],
        "--SVDPrefix",
        cfg["svd_prefix"],
        "--Output",
        output_prefix,
        "--NumPC",
        str(cfg.get("num_pc", 2)),
        "--NumThread",
        str(cfg.get("num_thread", 4)),
    ]
    if cfg.get("disable_sanity_check", False):
        cmd.append("--DisableSanityCheck")
    if cfg.get("output_pileup", False):
        cmd.append("--OutputPileup")
    known_af = cfg.get("known_af", "")
    if known_af:
        cmd.extend(["--KnownAF", known_af])
    extra_args = cfg.get("extra_args", [])
    if extra_args:
        cmd.extend(extra_args)
    return cmd


def normalize_contig(contig: str) -> str:
    contig = contig.strip()
    if contig.lower().startswith("chr"):
        contig = contig[3:]
    contig = contig.upper()
    if contig == "MT":
        contig = "M"
    return contig


def validate_panel(cfg: dict) -> None:
    svd_prefix = cfg["svd_prefix"]
    bed_path = Path(f"{svd_prefix}.bed")
    if not bed_path.exists():
        raise FileNotFoundError(f"Panel bed not found: {bed_path}")

    panel_build = cfg.get("panel_build", "")
    reference_build = cfg.get("reference_build", "")
    if not panel_build or not reference_build:
        raise ValueError("panel_build and reference_build must be set in config")
    if panel_build != reference_build:
        raise ValueError(
            f"Panel build ({panel_build}) does not match reference build ({reference_build})"
        )

    with bed_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            contig = normalize_contig(line.split("\t", 1)[0])
            if contig not in TARGET_CONTIGS:
                raise ValueError(
                    f"Non-target contig found in panel bed: {contig}"
                )


def parse_selfsm(selfsm_path: Path) -> dict:
    if not selfsm_path.exists():
        raise FileNotFoundError(f"SelfSM not found: {selfsm_path}")
    lines = [l.strip() for l in selfsm_path.read_text(encoding="utf-8").splitlines() if l.strip()]
    if len(lines) < 2:
        raise ValueError(f"SelfSM content is incomplete: {selfsm_path}")
    header = lines[0].lstrip("#").split("\t")
    values = lines[1].split("\t")
    if len(values) < len(header):
        raise ValueError(f"SelfSM row does not match header: {selfsm_path}")
    row = dict(zip(header, values))
    return {
        "freemix": float(row["FREEMIX"]),
        "freelk1": row.get("FREELK1", "NA"),
        "freelk0": row.get("FREELK0", "NA"),
        "chipmix": row.get("CHIPMIX", "NA"),
    }


def write_qc(output_prefix: Path, sample: str, metrics: dict) -> None:
    qc_path = Path(f"{output_prefix}.qc.tsv")
    status = "FAIL" if metrics["freemix"] > FREEMIX_FAIL_THRESHOLD else "PASS"
    header = ["sample", "freemix", "status", "freelk1", "freelk0", "chipmix"]
    row = [
        sample,
        f"{metrics['freemix']:.6f}",
        status,
        str(metrics["freelk1"]),
        str(metrics["freelk0"]),
        str(metrics["chipmix"]),
    ]
    qc_path.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")


def run(sample: str, bam: str, config: str, output_dir: str, log_dir: str, tmp_dir: str, force: bool) -> int:
    config_path = Path(config)
    cfg = load_config(config_path)
    validate_panel(cfg)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    log_path = Path(log_dir)
    log_path.mkdir(parents=True, exist_ok=True)
    sample_log = log_path / f"verifybamid2.{sample}.log"

    tmp_path = Path(tmp_dir)
    tmp_path.mkdir(parents=True, exist_ok=True)

    if log_has_success(sample_log) and not force:
        with sample_log.open("a", encoding="utf-8") as log_fh:
            log_fh.write(
                f"{dt.datetime.now().isoformat()}\tSTATUS=SUCCESS\tSKIP\n"
            )
        return 0

    bam_path = Path(bam)
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam}")

    output_prefix = str(output_path / sample)
    cmd = build_command(cfg, bam, output_prefix)

    with sample_log.open("a", encoding="utf-8") as log_fh:
        log_fh.write(f"{dt.datetime.now().isoformat()}\tSTART\n")
        log_fh.write(f"CMD: {' '.join(cmd)}\n")
        log_fh.flush()
        result = subprocess.run(cmd, stdout=log_fh, stderr=log_fh, check=False)
        if result.returncode == 0:
            log_fh.write(f"{dt.datetime.now().isoformat()}\tSTATUS=SUCCESS\n")
            metrics = parse_selfsm(Path(f"{output_prefix}.selfSM"))
            write_qc(Path(output_prefix), sample, metrics)
            log_fh.write(
                f"{dt.datetime.now().isoformat()}\tFREEMIX={metrics['freemix']:.6f}\n"
            )
            status = "FAIL" if metrics["freemix"] > FREEMIX_FAIL_THRESHOLD else "PASS"
            log_fh.write(
                f"{dt.datetime.now().isoformat()}\tQC_STATUS={status}\tTHRESHOLD={FREEMIX_FAIL_THRESHOLD:.2f}\n"
            )
        else:
            log_fh.write(
                f"{dt.datetime.now().isoformat()}\tSTATUS=FAIL\tCODE={result.returncode}\n"
            )
    return result.returncode


def main() -> None:
    parser = argparse.ArgumentParser(description="Run VerifyBamID2 for one sample.")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--bam", required=True, help="BAM/CRAM path")
    parser.add_argument("--config", required=True, help="Config YAML path")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--log-dir", required=True, help="Log directory")
    parser.add_argument("--tmp-dir", required=True, help="Temp directory")
    parser.add_argument("--force", action="store_true", help="Force rerun")
    args = parser.parse_args()

    raise SystemExit(
        run(
            sample=args.sample,
            bam=args.bam,
            config=args.config,
            output_dir=args.output_dir,
            log_dir=args.log_dir,
            tmp_dir=args.tmp_dir,
            force=args.force,
        )
    )


if __name__ == "__main__":
    main()
