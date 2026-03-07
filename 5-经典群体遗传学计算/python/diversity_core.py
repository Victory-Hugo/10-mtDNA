import inspect
from pathlib import Path
from typing import Iterable

import allel
import numpy as np
import pandas as pd
import zarr


WITHIN_ALLOWED = {"pi", "theta_w", "tajima_d", "haplotype_diversity"}
BETWEEN_ALLOWED = {"dxy", "hudson_fst"}


def parse_metrics(raw_value: str, allowed: set[str]) -> list[str]:
    metrics = [item.strip() for item in raw_value.split(",") if item.strip()]
    unknown = sorted(set(metrics) - allowed)
    if unknown:
        raise ValueError(f"Unsupported metrics: {', '.join(unknown)}")
    return metrics


def decode_array(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O"}:
        return np.array(
            [item.decode("utf-8") if isinstance(item, (bytes, bytearray)) else str(item) for item in values]
        )
    return values.astype(str)


def collapse_genotypes_to_haplotypes(gt: np.ndarray) -> tuple[np.ndarray, int]:
    if gt.ndim == 2:
        hap = gt.astype(np.int16, copy=False)
        return hap, 0
    if gt.ndim != 3:
        raise ValueError(f"Unsupported GT shape: {gt.shape}")
    if gt.shape[2] == 1:
        return gt[:, :, 0].astype(np.int16, copy=False), 0

    first = gt[:, :, 0].astype(np.int16, copy=False)
    second = gt[:, :, 1].astype(np.int16, copy=False)
    hap = first.copy()
    discordant = (first >= 0) & (second >= 0) & (first != second)
    missing = (first < 0) | (second < 0)
    hap[discordant | missing] = -1
    return hap, int(np.count_nonzero(discordant))


def invoke(func, **kwargs):
    signature = inspect.signature(func)
    accepted = {key: value for key, value in kwargs.items() if key in signature.parameters}
    return func(**accepted)


def safe_scalar(value) -> float:
    array = np.asarray(value)
    if array.size == 0:
        return float("nan")
    return float(array.reshape(-1)[0])


def build_windows(sequence_length: int, size: int, step: int) -> list[tuple[int, int]]:
    if size <= 0 or step <= 0:
        return []
    windows = []
    start = 1
    while start <= sequence_length:
        end = min(start + size, sequence_length + 1)
        windows.append((start, end))
        start += step
    return windows


def compute_haplotype_diversity(haplotype_array: allel.HaplotypeArray) -> float:
    if haplotype_array.shape[0] == 0 or haplotype_array.shape[1] == 0:
        return float("nan")
    try:
        return float(allel.haplotype_diversity(haplotype_array))
    except Exception:
        return float("nan")


def load_haplotype_dataset(zarr_path: str, contig_name: str, sequence_length: int) -> dict[str, object]:
    store_path = Path(zarr_path)
    if not store_path.exists():
        raise FileNotFoundError(f"Zarr store not found: {zarr_path}")

    root = zarr.open_group(str(store_path), mode="r")
    samples = decode_array(np.asarray(root["samples"][:]))
    chrom = decode_array(np.asarray(root["variants/CHROM"][:]))
    positions = np.asarray(root["variants/POS"][:], dtype=np.int64)
    gt = np.asarray(root["calldata/GT"][:])

    if contig_name:
        contig_mask = chrom == contig_name
    else:
        contig_mask = np.ones_like(positions, dtype=bool)
        contig_name = str(chrom[0]) if chrom.size else ""

    positions = positions[contig_mask]
    gt = gt[contig_mask]
    if sequence_length <= 0:
        sequence_length = int(positions.max()) if positions.size else 0

    hap, heterozygous_call_count = collapse_genotypes_to_haplotypes(gt)
    max_allele = int(hap[hap >= 0].max()) if np.any(hap >= 0) else 1

    return {
        "samples": samples,
        "positions": positions,
        "hap": hap,
        "max_allele": max_allele,
        "contig_name": contig_name,
        "sequence_length": sequence_length,
        "heterozygous_call_count": heterozygous_call_count,
    }


def build_group_indices(
    metadata: pd.DataFrame,
    samples: np.ndarray,
    sample_id_column: str,
    group_column: str,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    ordered_metadata = metadata.set_index(sample_id_column).loc[samples].reset_index()
    sample_to_index = {sample_id: idx for idx, sample_id in enumerate(samples.tolist())}
    group_to_indices = {
        group: np.array([sample_to_index[sample_id] for sample_id in frame[sample_id_column]], dtype=np.int64)
        for group, frame in ordered_metadata.groupby(group_column, sort=True)
    }
    return ordered_metadata, group_to_indices


def compute_within_group_data(
    indices: np.ndarray,
    hap: np.ndarray,
    positions: np.ndarray,
    sequence_length: int,
    max_allele: int,
    metrics: Iterable[str],
) -> tuple[dict[str, float], dict[str, float], allel.AlleleCountsArray, allel.HaplotypeArray]:
    group_hap = allel.HaplotypeArray(hap[:, indices])
    allele_counts = group_hap.count_alleles(max_allele=max_allele)
    segregating = allele_counts.is_segregating()
    callable_sites = np.any(group_hap >= 0, axis=1)

    info = {
        "n_samples": int(len(indices)),
        "n_sites_used": int(sequence_length),
        "n_variant_sites": int(np.count_nonzero(callable_sites)),
        "n_segregating_sites": int(np.count_nonzero(segregating)),
        "missing_rate_mean": float(np.mean(group_hap < 0)),
    }

    metric_values: dict[str, float] = {}
    if "pi" in metrics:
        metric_values["pi"] = safe_scalar(
            invoke(allel.sequence_diversity, pos=positions, ac=allele_counts, start=1, stop=sequence_length + 1)
        )
    if "theta_w" in metrics:
        metric_values["theta_w"] = safe_scalar(
            invoke(allel.watterson_theta, pos=positions, ac=allele_counts, start=1, stop=sequence_length + 1)
        )
    if "tajima_d" in metrics:
        try:
            metric_values["tajima_d"] = safe_scalar(
                invoke(allel.tajima_d, ac=allele_counts, pos=positions, start=1, stop=sequence_length + 1)
            )
        except Exception:
            metric_values["tajima_d"] = float("nan")
    if "haplotype_diversity" in metrics:
        metric_values["haplotype_diversity"] = compute_haplotype_diversity(group_hap)

    return info, metric_values, allele_counts, group_hap


def compute_between_group_data(
    ac1: allel.AlleleCountsArray,
    ac2: allel.AlleleCountsArray,
    positions: np.ndarray,
    sequence_length: int,
    metrics: Iterable[str],
) -> tuple[dict[str, float], dict[str, float]]:
    called_mask = (ac1.sum(axis=1) > 0) & (ac2.sum(axis=1) > 0)
    info = {
        "n_sites_used": int(sequence_length),
        "shared_variant_sites": int(np.count_nonzero(called_mask)),
    }

    metric_values: dict[str, float] = {}
    if "dxy" in metrics:
        metric_values["dxy"] = safe_scalar(
            invoke(allel.sequence_divergence, pos=positions, ac1=ac1, ac2=ac2, start=1, stop=sequence_length + 1)
        )
    if "hudson_fst" in metrics:
        numerator, denominator = allel.hudson_fst(ac1, ac2)
        denominator_sum = np.nansum(denominator)
        metric_values["hudson_fst"] = float(np.nansum(numerator) / denominator_sum) if denominator_sum > 0 else float(
            "nan"
        )

    return info, metric_values
