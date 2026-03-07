import argparse
import logging
from itertools import combinations
from pathlib import Path
from typing import Iterable

import allel
import numpy as np
import pandas as pd

from diversity_core import (
    BETWEEN_ALLOWED,
    WITHIN_ALLOWED,
    build_group_indices,
    build_windows,
    compute_between_group_data,
    compute_haplotype_diversity,
    compute_within_group_data,
    invoke,
    load_haplotype_dataset,
    parse_metrics,
    safe_scalar,
)


log = logging.getLogger(__name__)


def compute_within_windows(
    group: str,
    group_hap: allel.HaplotypeArray,
    allele_counts: allel.AlleleCountsArray,
    positions: np.ndarray,
    windows: list[tuple[int, int]],
    metrics: Iterable[str],
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for start, end in windows:
        mask = (positions >= start) & (positions < end)
        idx = np.where(mask)[0]
        pos_window = positions[idx]
        ac_window = allele_counts[idx]
        hap_window = group_hap[idx]
        n_variants = int(idx.size)
        accessible_bases = int(end - start)
        metric_values: dict[str, float] = {}
        if "pi" in metrics:
            metric_values["pi"] = float("nan") if n_variants == 0 else safe_scalar(
                invoke(allel.sequence_diversity, pos=pos_window, ac=ac_window, start=start, stop=end)
            )
        if "theta_w" in metrics:
            metric_values["theta_w"] = float("nan") if n_variants == 0 else safe_scalar(
                invoke(allel.watterson_theta, pos=pos_window, ac=ac_window, start=start, stop=end)
            )
        if "tajima_d" in metrics:
            if n_variants == 0:
                metric_values["tajima_d"] = float("nan")
            else:
                try:
                    metric_values["tajima_d"] = safe_scalar(
                        invoke(allel.tajima_d, ac=ac_window, pos=pos_window, start=start, stop=end)
                    )
                except Exception:
                    metric_values["tajima_d"] = float("nan")
        if "haplotype_diversity" in metrics:
            metric_values["haplotype_diversity"] = compute_haplotype_diversity(hap_window)

        for metric_name, metric_value in metric_values.items():
            rows.append(
                {
                    "contig": "",
                    "start": start,
                    "end": end - 1,
                    "group_or_pair": group,
                    "metric": metric_name,
                    "value": metric_value,
                    "n_variants": n_variants,
                    "n_accessible_bases": accessible_bases,
                }
            )
    return rows


def compute_between_summary(
    group1: str,
    group2: str,
    ac1: allel.AlleleCountsArray,
    ac2: allel.AlleleCountsArray,
    positions: np.ndarray,
    sequence_length: int,
    metrics: Iterable[str],
) -> dict[str, float]:
    info, metric_values = compute_between_group_data(
        ac1=ac1,
        ac2=ac2,
        positions=positions,
        sequence_length=sequence_length,
        metrics=metrics,
    )
    row: dict[str, float] = {"group1": group1, "group2": group2, **info}
    row.update(metric_values)
    return row


def compute_between_windows(
    group1: str,
    group2: str,
    ac1: allel.AlleleCountsArray,
    ac2: allel.AlleleCountsArray,
    positions: np.ndarray,
    windows: list[tuple[int, int]],
    metrics: Iterable[str],
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    pair_name = f"{group1}|{group2}"
    for start, end in windows:
        mask = (positions >= start) & (positions < end)
        idx = np.where(mask)[0]
        pos_window = positions[idx]
        ac1_window = ac1[idx]
        ac2_window = ac2[idx]
        n_variants = int(idx.size)
        accessible_bases = int(end - start)
        metric_values: dict[str, float] = {}
        if "dxy" in metrics:
            metric_values["dxy"] = float("nan") if n_variants == 0 else safe_scalar(
                invoke(allel.sequence_divergence, pos=pos_window, ac1=ac1_window, ac2=ac2_window, start=start, stop=end)
            )
        if "hudson_fst" in metrics:
            if n_variants == 0:
                metric_values["hudson_fst"] = float("nan")
            else:
                numerator, denominator = allel.hudson_fst(ac1_window, ac2_window)
                denominator_sum = np.nansum(denominator)
                metric_values["hudson_fst"] = (
                    float(np.nansum(numerator) / denominator_sum) if denominator_sum > 0 else float("nan")
                )

        for metric_name, metric_value in metric_values.items():
            rows.append(
                {
                    "contig": "",
                    "start": start,
                    "end": end - 1,
                    "group_or_pair": pair_name,
                    "metric": metric_name,
                    "value": metric_value,
                    "n_variants": n_variants,
                    "n_accessible_bases": accessible_bases,
                }
            )
    return rows


def run(
    zarr_path: str,
    matched_samples: str,
    sample_id_column: str,
    group_column: str,
    contig_name: str,
    sequence_length: int,
    within_metrics: str,
    between_metrics: str,
    window_enable: int,
    window_size: int,
    window_step: int,
    within_output: str,
    between_output: str,
    window_output: str,
) -> dict[str, int]:
    within_metric_list = parse_metrics(within_metrics, WITHIN_ALLOWED)
    between_metric_list = parse_metrics(between_metrics, BETWEEN_ALLOWED)
    metadata = pd.read_csv(matched_samples, sep="\t", dtype=str)
    dataset = load_haplotype_dataset(zarr_path=zarr_path, contig_name=contig_name, sequence_length=sequence_length)
    samples = dataset["samples"]
    positions = dataset["positions"]
    hap = dataset["hap"]
    max_allele = dataset["max_allele"]
    contig_name = dataset["contig_name"]
    sequence_length = dataset["sequence_length"]
    heterozygous_call_count = dataset["heterozygous_call_count"]

    windows = build_windows(sequence_length, window_size, window_step) if window_enable else []
    within_rows: list[dict[str, float]] = []
    between_rows: list[dict[str, float]] = []
    window_rows: list[dict[str, float]] = []
    group_allele_counts: dict[str, allel.AlleleCountsArray] = {}

    Path(within_output).parent.mkdir(parents=True, exist_ok=True)
    Path(between_output).parent.mkdir(parents=True, exist_ok=True)
    Path(window_output).parent.mkdir(parents=True, exist_ok=True)

    _, group_to_indices = build_group_indices(
        metadata=metadata,
        samples=samples,
        sample_id_column=sample_id_column,
        group_column=group_column,
    )

    for index, (group, indices) in enumerate(group_to_indices.items(), start=1):
        info, metric_values, allele_counts, group_hap = compute_within_group_data(
            indices=indices,
            hap=hap,
            positions=positions,
            sequence_length=sequence_length,
            max_allele=max_allele,
            metrics=within_metric_list,
        )
        row = {
            "group": group,
            **info,
            "heterozygous_calls_masked_global": heterozygous_call_count,
            **metric_values,
        }
        within_rows.append(row)
        group_allele_counts[group] = allele_counts
        if windows:
            window_rows.extend(compute_within_windows(group, group_hap, allele_counts, positions, windows, within_metric_list))
        if index % 25 == 0 or index == len(group_to_indices):
            log.info("Computed within-group metrics for %d/%d groups", index, len(group_to_indices))

    group_names = list(group_to_indices.keys())
    pair_total = len(group_names) * (len(group_names) - 1) // 2
    for pair_index, (group1, group2) in enumerate(combinations(group_names, 2), start=1):
        ac1 = group_allele_counts[group1]
        ac2 = group_allele_counts[group2]
        between_rows.append(
            compute_between_summary(group1, group2, ac1, ac2, positions, sequence_length, between_metric_list)
        )
        if windows:
            window_rows.extend(compute_between_windows(group1, group2, ac1, ac2, positions, windows, between_metric_list))
        if pair_index % 500 == 0 or pair_index == pair_total:
            log.info("Computed between-group metrics for %d/%d pairs", pair_index, pair_total)

    within_df = pd.DataFrame(within_rows).sort_values("group").reset_index(drop=True)
    between_df = pd.DataFrame(between_rows).sort_values(["group1", "group2"]).reset_index(drop=True)
    within_df.to_csv(within_output, sep="\t", index=False)
    between_df.to_csv(between_output, sep="\t", index=False)

    if windows:
        window_df = pd.DataFrame(window_rows)
        window_df["contig"] = contig_name
        window_df.to_csv(window_output, sep="\t", index=False)
    else:
        pd.DataFrame(
            columns=["contig", "start", "end", "group_or_pair", "metric", "value", "n_variants", "n_accessible_bases"]
        ).to_csv(window_output, sep="\t", index=False)

    log.info("Wrote %d within-group rows and %d between-group rows", len(within_df), len(between_df))
    return {
        "within_groups": len(within_df),
        "between_pairs": len(between_df),
        "window_rows": len(window_rows),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute within-group and between-group diversity metrics from Zarr.")
    parser.add_argument("--zarr-path", required=True)
    parser.add_argument("--matched-samples", required=True)
    parser.add_argument("--sample-id-column", required=True)
    parser.add_argument("--group-column", required=True)
    parser.add_argument("--contig-name", required=True)
    parser.add_argument("--sequence-length", type=int, default=0)
    parser.add_argument("--within-metrics", required=True)
    parser.add_argument("--between-metrics", required=True)
    parser.add_argument("--window-enable", type=int, default=0)
    parser.add_argument("--window-size", type=int, default=0)
    parser.add_argument("--window-step", type=int, default=0)
    parser.add_argument("--within-output", required=True)
    parser.add_argument("--between-output", required=True)
    parser.add_argument("--window-output", required=True)
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
        window_enable=args.window_enable,
        window_size=args.window_size,
        window_step=args.window_step,
        within_output=args.within_output,
        between_output=args.between_output,
        window_output=args.window_output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
