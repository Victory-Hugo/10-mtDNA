import argparse
import json
import logging
import shlex
from pathlib import Path


log = logging.getLogger(__name__)


def resolve_path(project_root: Path, value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        return str(path)
    return str((project_root / path).resolve())


def load_config(config_path: str, project_root: str | None = None) -> dict[str, str]:
    config_file = Path(config_path).resolve()
    root = Path(project_root).resolve() if project_root else config_file.parent.parent.resolve()
    with config_file.open("r", encoding="utf-8") as handle:
        config = json.load(handle)

    metrics = config.get("metrics", {})
    windows = config.get("windows", {})
    runtime = config.get("runtime", {})
    resampling = config.get("resampling", {})

    project_name = config.get("project_name", "scikit_allel_pipeline")
    output_dir = resolve_path(root, config.get("output_dir", "output"))
    tmp_dir = resolve_path(root, config.get("tmp_dir", "data/tmp"))
    meta_dir = resolve_path(root, config.get("meta_dir", "meta"))
    log_dir = resolve_path(root, config.get("log_dir", "logs"))
    python_bin = resolve_path(root, runtime.get("python_bin", "python3"))

    env = {
        "PROJECT_ROOT": str(root),
        "PROJECT_NAME": project_name,
        "INPUT_VCF": resolve_path(root, config["input_vcf"]),
        "SAMPLE_TABLE": resolve_path(root, config["sample_table"]),
        "SAMPLE_ID_COLUMN": config.get("sample_id_column", "ID"),
        "GROUP_COLUMN": config.get("group_column", "Group"),
        "MIN_GROUP_SIZE": str(config.get("min_group_size", 1)),
        "CONTIG_NAME": config.get("contig_name", ""),
        "SEQUENCE_LENGTH": str(config.get("sequence_length", 0)),
        "OUTPUT_DIR": output_dir,
        "TMP_DIR": tmp_dir,
        "META_DIR": meta_dir,
        "LOG_DIR": log_dir,
        "MATCHED_SAMPLES_TSV": str(Path(meta_dir) / "matched_samples.tsv"),
        "GROUP_COUNTS_TSV": str(Path(meta_dir) / "group_counts.tsv"),
        "INPUT_SUMMARY_TSV": str(Path(output_dir) / "0-qc" / "input_summary.tsv"),
        "WITHIN_SUMMARY_TSV": str(Path(output_dir) / "1-within" / "within_group_summary.tsv"),
        "BETWEEN_SUMMARY_TSV": str(Path(output_dir) / "2-between" / "between_group_summary.tsv"),
        "WINDOW_METRICS_TSV": str(Path(output_dir) / "3-windows" / "windowed_metrics.tsv"),
        "BOOTSTRAP_DIR": str(Path(output_dir) / "4-bootstrap"),
        "WITHIN_BOOTSTRAP_SUMMARY_TSV": str(Path(output_dir) / "4-bootstrap" / "within_group_bootstrap_summary.tsv"),
        "BETWEEN_BOOTSTRAP_SUMMARY_TSV": str(Path(output_dir) / "4-bootstrap" / "between_group_bootstrap_summary.tsv"),
        "BOOTSTRAP_RUN_SUMMARY_TSV": str(Path(output_dir) / "4-bootstrap" / "bootstrap_run_summary.tsv"),
        "WRITE_BOOTSTRAP_REPLICATE_TABLES": "1" if resampling.get("write_replicate_tables", False) else "0",
        "WITHIN_BOOTSTRAP_REPLICATES_TSV": str(
            Path(output_dir) / "4-bootstrap" / "within_group_bootstrap_replicates.tsv"
        ),
        "BETWEEN_BOOTSTRAP_REPLICATES_TSV": str(
            Path(output_dir) / "4-bootstrap" / "between_group_bootstrap_replicates.tsv"
        ),
        "ZARR_PATH": str(Path(tmp_dir) / f"{project_name}.zarr"),
        "WITHIN_METRICS": ",".join(metrics.get("within_groups", [])),
        "BETWEEN_METRICS": ",".join(metrics.get("between_groups", [])),
        "WINDOW_ENABLE": "1" if windows.get("enable", False) else "0",
        "WINDOW_SIZE": str(windows.get("size", 0)),
        "WINDOW_STEP": str(windows.get("step", 0)),
        "RESAMPLING_ENABLE": "1" if resampling.get("enable", False) else "0",
        "RESAMPLE_STRATEGY": resampling.get("sample_size_strategy", "min_group_size"),
        "RESAMPLING_DOWNSAMPLE": str(resampling.get("resampling_downsample", 20)),
        "BOOTSTRAP_REPLICATES": str(resampling.get("bootstrap_replicates", 1000)),
        "BOOTSTRAP_RANDOM_SEED": str(resampling.get("random_seed", 42)),
        "PYTHON_BIN": python_bin,
        "N_THREADS": str(runtime.get("n_threads", 1)),
        "OVERWRITE_ZARR": "1" if runtime.get("overwrite_zarr", False) else "0",
        "ZARR_CHUNK_LENGTH": str(runtime.get("zarr_chunk_length", 2048)),
        "ZARR_CHUNK_WIDTH": str(runtime.get("zarr_chunk_width", 128)),
        "ALT_NUMBER": str(runtime.get("alt_number", 4)),
    }
    return env


def run(config: str, project_root: str | None = None) -> dict[str, str]:
    return load_config(config_path=config, project_root=project_root)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Read Config.json and emit shell variables.")
    parser.add_argument("--config", required=True, help="Path to conf/Config.json")
    parser.add_argument("--project-root", help="Optional project root override")
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    env = run(config=args.config, project_root=args.project_root)
    for key, value in env.items():
        print(f"{key}={shlex.quote(str(value))}")
    log.info("Emitted %d config variables", len(env))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
