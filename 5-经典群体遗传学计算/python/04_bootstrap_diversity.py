import argparse
import logging
import math
import multiprocessing as mp
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from diversity_core import (
    BETWEEN_ALLOWED,
    WITHIN_ALLOWED,
    build_group_indices,
    compute_between_group_data,
    compute_within_group_data,
    load_haplotype_dataset,
    parse_metrics,
)


log = logging.getLogger(__name__)
BOOTSTRAP_CONTEXT: dict[str, object] = {}


def _set_bootstrap_context(context: dict[str, object]) -> None:
    global BOOTSTRAP_CONTEXT
    BOOTSTRAP_CONTEXT = context


def _run_bootstrap_chunk(replicate_ids: list[int]) -> tuple[list[int], np.ndarray, np.ndarray]:
    hap = BOOTSTRAP_CONTEXT["hap"]
    positions = BOOTSTRAP_CONTEXT["positions"]
    sequence_length = BOOTSTRAP_CONTEXT["sequence_length"]
    max_allele = BOOTSTRAP_CONTEXT["max_allele"]
    group_names = BOOTSTRAP_CONTEXT["group_names"]
    group_to_indices = BOOTSTRAP_CONTEXT["group_to_indices"]
    within_metrics = BOOTSTRAP_CONTEXT["within_metrics"]
    between_metrics = BOOTSTRAP_CONTEXT["between_metrics"]
    group_pairs = BOOTSTRAP_CONTEXT["group_pairs"]
    resample_n = BOOTSTRAP_CONTEXT["resample_n"]
    random_seed = BOOTSTRAP_CONTEXT["random_seed"]

    within_chunk = np.full((len(group_names), len(within_metrics), len(replicate_ids)), np.nan, dtype=np.float32)
    between_chunk = np.full((len(group_pairs), len(between_metrics), len(replicate_ids)), np.nan, dtype=np.float32)

    for local_rep_index, replicate_id in enumerate(replicate_ids):
        rng = np.random.default_rng(random_seed + replicate_id)
        bootstrap_acs = {}
        for group_index, group_name in enumerate(group_names):
            sampled_indices = rng.choice(group_to_indices[group_name], size=resample_n, replace=True)
            _, metric_values, allele_counts, _ = compute_within_group_data(
                indices=sampled_indices,
                hap=hap,
                positions=positions,
                sequence_length=sequence_length,
                max_allele=max_allele,
                metrics=within_metrics,
            )
            bootstrap_acs[group_name] = allele_counts
            for metric_index, metric_name in enumerate(within_metrics):
                within_chunk[group_index, metric_index, local_rep_index] = metric_values.get(metric_name, np.nan)

        for pair_index, (group1, group2) in enumerate(group_pairs):
            _, metric_values = compute_between_group_data(
                bootstrap_acs[group1],
                bootstrap_acs[group2],
                positions=positions,
                sequence_length=sequence_length,
                metrics=between_metrics,
            )
            for metric_index, metric_name in enumerate(between_metrics):
                between_chunk[pair_index, metric_index, local_rep_index] = metric_values.get(metric_name, np.nan)

    return replicate_ids, within_chunk, between_chunk


def choose_resample_size(group_to_indices: dict[str, np.ndarray], strategy: str) -> int:
    if strategy != "min_group_size":
        raise ValueError(f"Unsupported sample_size_strategy: {strategy}")
    if len(group_to_indices) < 2:
        raise ValueError("At least two groups are required for bootstrap analysis.")
    return min(int(len(indices)) for indices in group_to_indices.values())


def load_reference_within(path: str, metrics: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    missing_metrics = sorted(set(metrics) - set(frame.columns))
    if missing_metrics:
        raise ValueError(f"Within reference missing metrics: {', '.join(missing_metrics)}")
    melted_frames = []
    for metric_name in metrics:
        metric_frame = frame[["group", "n_samples", metric_name]].copy()
        metric_frame["metric"] = metric_name
        metric_frame = metric_frame.rename(columns={"n_samples": "original_n_samples", metric_name: "full_data_value"})
        melted_frames.append(metric_frame)
    return pd.concat(melted_frames, ignore_index=True)


def load_reference_between(path: str, metrics: list[str], group_size_map: dict[str, int]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    missing_metrics = sorted(set(metrics) - set(frame.columns))
    if missing_metrics:
        raise ValueError(f"Between reference missing metrics: {', '.join(missing_metrics)}")
    melted_frames = []
    for metric_name in metrics:
        metric_frame = frame[["group1", "group2", metric_name]].copy()
        metric_frame["metric"] = metric_name
        metric_frame["original_n_samples_group1"] = metric_frame["group1"].map(group_size_map)
        metric_frame["original_n_samples_group2"] = metric_frame["group2"].map(group_size_map)
        metric_frame = metric_frame.rename(columns={metric_name: "full_data_value"})
        melted_frames.append(metric_frame)
    output = pd.concat(melted_frames, ignore_index=True)
    if output[["original_n_samples_group1", "original_n_samples_group2"]].isna().any().any():
        raise ValueError("Unable to map original group sizes for between-group bootstrap output.")
    return output


def summarize_bootstrap_array(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with np.errstate(invalid="ignore"):
        valid_counts = np.sum(~np.isnan(values), axis=2)
        mean_values = np.nanmean(values, axis=2, dtype=np.float64)
        std_values = np.nanstd(values, axis=2, dtype=np.float64)
        p2_5 = np.nanpercentile(values, 2.5, axis=2)
        p50 = np.nanpercentile(values, 50.0, axis=2)
        p97_5 = np.nanpercentile(values, 97.5, axis=2)
    return valid_counts, mean_values, std_values, p2_5, p50, p97_5


def build_chunks(total_replicates: int, n_threads: int) -> list[list[int]]:
    if total_replicates < 1:
        return []
    n_chunks = min(max(n_threads, 1), total_replicates)
    chunk_size = math.ceil(total_replicates / n_chunks)
    return [list(range(start, min(start + chunk_size, total_replicates))) for start in range(0, total_replicates, chunk_size)]


def run(
    zarr_path: str,
    matched_samples: str,
    sample_id_column: str,
    group_column: str,
    contig_name: str,
    sequence_length: int,
    within_metrics: str,
    between_metrics: str,
    bootstrap_replicates: int,
    random_seed: int,
    sample_size_strategy: str,
    within_reference: str,
    between_reference: str,
    within_bootstrap_output: str,
    between_bootstrap_output: str,
    run_summary_output: str,
    n_threads: int,
) -> dict[str, int]:
    if bootstrap_replicates < 1:
        raise ValueError("bootstrap_replicates must be >= 1")

    within_metric_list = parse_metrics(within_metrics, WITHIN_ALLOWED)
    between_metric_list = parse_metrics(between_metrics, BETWEEN_ALLOWED)
    metadata = pd.read_csv(matched_samples, sep="\t", dtype=str)
    dataset = load_haplotype_dataset(zarr_path=zarr_path, contig_name=contig_name, sequence_length=sequence_length)
    _, group_to_indices = build_group_indices(
        metadata=metadata,
        samples=dataset["samples"],
        sample_id_column=sample_id_column,
        group_column=group_column,
    )

    resample_n = choose_resample_size(group_to_indices, sample_size_strategy)
    group_names = sorted(group_to_indices)
    group_pairs = list(combinations(group_names, 2))
    group_size_map = {group_name: int(len(group_to_indices[group_name])) for group_name in group_names}

    Path(within_bootstrap_output).parent.mkdir(parents=True, exist_ok=True)
    Path(between_bootstrap_output).parent.mkdir(parents=True, exist_ok=True)
    Path(run_summary_output).parent.mkdir(parents=True, exist_ok=True)

    within_values = np.full((len(group_names), len(within_metric_list), bootstrap_replicates), np.nan, dtype=np.float32)
    between_values = np.full(
        (len(group_pairs), len(between_metric_list), bootstrap_replicates),
        np.nan,
        dtype=np.float32,
    )

    chunks = build_chunks(bootstrap_replicates, max(1, n_threads))
    bootstrap_context = {
        "hap": dataset["hap"],
        "positions": dataset["positions"],
        "sequence_length": dataset["sequence_length"],
        "max_allele": dataset["max_allele"],
        "group_names": group_names,
        "group_to_indices": group_to_indices,
        "within_metrics": within_metric_list,
        "between_metrics": between_metric_list,
        "group_pairs": group_pairs,
        "resample_n": resample_n,
        "random_seed": int(random_seed),
    }

    use_parallel = max(1, n_threads) > 1 and len(chunks) > 1 and "fork" in mp.get_all_start_methods()
    if use_parallel:
        _set_bootstrap_context(bootstrap_context)
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=min(max(1, n_threads), len(chunks))) as pool:
            for completed_index, (replicate_ids, within_chunk, between_chunk) in enumerate(
                pool.imap_unordered(_run_bootstrap_chunk, chunks),
                start=1,
            ):
                within_values[:, :, replicate_ids] = within_chunk
                between_values[:, :, replicate_ids] = between_chunk
                log.info(
                    "Completed bootstrap chunk %d/%d (%d replicates)",
                    completed_index,
                    len(chunks),
                    len(replicate_ids),
                )
    else:
        _set_bootstrap_context(bootstrap_context)
        for completed_index, chunk in enumerate(chunks, start=1):
            replicate_ids, within_chunk, between_chunk = _run_bootstrap_chunk(chunk)
            within_values[:, :, replicate_ids] = within_chunk
            between_values[:, :, replicate_ids] = between_chunk
            log.info(
                "Completed bootstrap chunk %d/%d (%d replicates)",
                completed_index,
                len(chunks),
                len(replicate_ids),
            )

    within_reference_frame = load_reference_within(within_reference, within_metric_list)
    between_reference_frame = load_reference_between(between_reference, between_metric_list, group_size_map)

    within_valid, within_mean, within_std, within_p2_5, within_p50, within_p97_5 = summarize_bootstrap_array(within_values)
    between_valid, between_mean, between_std, between_p2_5, between_p50, between_p97_5 = summarize_bootstrap_array(
        between_values
    )

    within_reference_map = {(row.group, row.metric): row for row in within_reference_frame.itertuples(index=False)}
    between_reference_map = {
        (row.group1, row.group2, row.metric): row for row in between_reference_frame.itertuples(index=False)
    }

    within_rows = []
    for group_index, group_name in enumerate(group_names):
        for metric_index, metric_name in enumerate(within_metric_list):
            reference_row = within_reference_map[(group_name, metric_name)]
            within_rows.append(
                {
                    "group": group_name,
                    "original_n_samples": int(reference_row.original_n_samples),
                    "resample_n": resample_n,
                    "bootstrap_replicates_requested": bootstrap_replicates,
                    "n_valid_bootstrap": int(within_valid[group_index, metric_index]),
                    "metric": metric_name,
                    "full_data_value": float(reference_row.full_data_value),
                    "bootstrap_mean": float(within_mean[group_index, metric_index]),
                    "bootstrap_sd": float(within_std[group_index, metric_index]),
                    "bootstrap_p2_5": float(within_p2_5[group_index, metric_index]),
                    "bootstrap_p50": float(within_p50[group_index, metric_index]),
                    "bootstrap_p97_5": float(within_p97_5[group_index, metric_index]),
                }
            )

    between_rows = []
    for pair_index, (group1, group2) in enumerate(group_pairs):
        for metric_index, metric_name in enumerate(between_metric_list):
            reference_row = between_reference_map[(group1, group2, metric_name)]
            between_rows.append(
                {
                    "group1": group1,
                    "group2": group2,
                    "original_n_samples_group1": int(reference_row.original_n_samples_group1),
                    "original_n_samples_group2": int(reference_row.original_n_samples_group2),
                    "resample_n": resample_n,
                    "bootstrap_replicates_requested": bootstrap_replicates,
                    "n_valid_bootstrap": int(between_valid[pair_index, metric_index]),
                    "metric": metric_name,
                    "full_data_value": float(reference_row.full_data_value),
                    "bootstrap_mean": float(between_mean[pair_index, metric_index]),
                    "bootstrap_sd": float(between_std[pair_index, metric_index]),
                    "bootstrap_p2_5": float(between_p2_5[pair_index, metric_index]),
                    "bootstrap_p50": float(between_p50[pair_index, metric_index]),
                    "bootstrap_p97_5": float(between_p97_5[pair_index, metric_index]),
                }
            )

    pd.DataFrame(within_rows).to_csv(within_bootstrap_output, sep="\t", index=False)
    pd.DataFrame(between_rows).to_csv(between_bootstrap_output, sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "n_groups": len(group_names),
                "n_pairs": len(group_pairs),
                "resample_n": resample_n,
                "sample_size_strategy": sample_size_strategy,
                "bootstrap_enabled": 1,
                "bootstrap_replicates_requested": bootstrap_replicates,
                "random_seed": random_seed,
                "within_metric_count": len(within_metric_list),
                "between_metric_count": len(between_metric_list),
                "window_bootstrap_enabled": 0,
            }
        ]
    ).to_csv(run_summary_output, sep="\t", index=False)

    log.info(
        "Wrote bootstrap summaries for %d groups and %d group pairs using resample_n=%d",
        len(group_names),
        len(group_pairs),
        resample_n,
    )
    return {
        "within_rows": len(within_rows),
        "between_rows": len(between_rows),
        "resample_n": resample_n,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run group-wise sample bootstrap on diversity metrics.")
    parser.add_argument("--zarr-path", required=True)
    parser.add_argument("--matched-samples", required=True)
    parser.add_argument("--sample-id-column", required=True)
    parser.add_argument("--group-column", required=True)
    parser.add_argument("--contig-name", required=True)
    parser.add_argument("--sequence-length", type=int, default=0)
    parser.add_argument("--within-metrics", required=True)
    parser.add_argument("--between-metrics", required=True)
    parser.add_argument("--bootstrap-replicates", type=int, required=True)
    parser.add_argument("--random-seed", type=int, required=True)
    parser.add_argument("--sample-size-strategy", required=True)
    parser.add_argument("--within-reference", required=True)
    parser.add_argument("--between-reference", required=True)
    parser.add_argument("--within-bootstrap-output", required=True)
    parser.add_argument("--between-bootstrap-output", required=True)
    parser.add_argument("--run-summary-output", required=True)
    parser.add_argument("--n-threads", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = build_parser().parse_args(argv)
    run(
        zarr_path=args.zarr_path,
        matched_samples=args.matched_samples,
        sample_id_column=args.sample_id_column,
        group_column=args.group_column,
        contig_name=args.contig_name,
        sequence_length=args.sequence_length,
        within_metrics=args.within_metrics,
        between_metrics=args.between_metrics,
        bootstrap_replicates=args.bootstrap_replicates,
        random_seed=args.random_seed,
        sample_size_strategy=args.sample_size_strategy,
        within_reference=args.within_reference,
        between_reference=args.between_reference,
        within_bootstrap_output=args.within_bootstrap_output,
        between_bootstrap_output=args.between_bootstrap_output,
        run_summary_output=args.run_summary_output,
        n_threads=args.n_threads,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
