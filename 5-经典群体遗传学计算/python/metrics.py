from __future__ import annotations

from typing import Dict, List, Tuple, TypedDict
import logging

import numpy as np
import pandas as pd
import allel

logger = logging.getLogger(__name__)


class SfsResult(TypedDict):
    sfs: np.ndarray | None
    n_samples: int
    n_segregating: int


def calculate_pi(gt_array: allel.GenotypeArray, populations: Dict[str, List[int]], pos: np.ndarray) -> Dict[str, float]:
    pi_values = {}
    logger.info("计算 π: 群体数 %s", len(populations))
    for pop_name, pop_indices in populations.items():
        pop_gt = gt_array.take(pop_indices, axis=1)
        ac = pop_gt.count_alleles()
        pi_values[pop_name] = allel.sequence_diversity(pos, ac)
    return pi_values


def calculate_theta_w(gt_array: allel.GenotypeArray, populations: Dict[str, List[int]], pos: np.ndarray) -> Dict[str, float]:
    theta_w_values = {}
    logger.info("计算 θw: 群体数 %s", len(populations))
    for pop_name, pop_indices in populations.items():
        pop_gt = gt_array.take(pop_indices, axis=1)
        ac = pop_gt.count_alleles()
        theta_w_values[pop_name] = allel.watterson_theta(pos, ac)
    return theta_w_values


def calculate_tajima_d(gt_array: allel.GenotypeArray, populations: Dict[str, List[int]], pos: np.ndarray) -> Dict[str, float]:
    tajima_d_values = {}
    logger.info("计算 Tajima's D: 群体数 %s", len(populations))
    for pop_name, pop_indices in populations.items():
        pop_gt = gt_array.take(pop_indices, axis=1)
        ac = pop_gt.count_alleles()
        tajima_d_values[pop_name] = allel.tajima_d(ac, pos=pos)
    return tajima_d_values


def calculate_sfs(
    gt_array: allel.GenotypeArray,
    populations: Dict[str, List[int]],
    folded: bool = True,
) -> Dict[str, SfsResult]:
    sfs_values = {}
    logger.info("计算 SFS: 群体数 %s, folded=%s", len(populations), folded)
    for pop_name, pop_indices in populations.items():
        pop_gt = gt_array.take(pop_indices, axis=1)
        ac = pop_gt.count_alleles()
        n_samples = pop_gt.n_samples

        if ac.shape[1] > 2:
            logger.warning("群体 %s 存在多等位基因位点，SFS 将仅使用二等位位点", pop_name)

        if hasattr(ac, "is_biallelic"):
            biallelic = ac.is_biallelic()
        else:
            biallelic = ac.max_allele() <= 1

        if hasattr(ac, "is_segregating"):
            segregating = ac.is_segregating()
        else:
            segregating = np.sum(ac > 0, axis=1) > 1

        ac_seg = ac[biallelic & segregating]
        n_segregating = int(ac_seg.shape[0])

        if len(ac_seg) == 0:
            logger.warning("群体 %s 没有分离位点", pop_name)
            sfs_values[pop_name] = {
                "sfs": None,
                "n_samples": int(n_samples),
                "n_segregating": n_segregating,
            }
            continue

        try:
            if folded:
                sfs_values[pop_name] = {
                    "sfs": allel.sfs_folded(ac_seg),
                    "n_samples": int(n_samples),
                    "n_segregating": n_segregating,
                }
            else:
                dac = ac_seg[:, 1]
                sfs_values[pop_name] = {
                    "sfs": allel.sfs(dac),
                    "n_samples": int(n_samples),
                    "n_segregating": n_segregating,
                }
        except ValueError as exc:
            logger.warning("群体 %s SFS 计算失败: %s", pop_name, exc)
            sfs_values[pop_name] = {
                "sfs": None,
                "n_samples": int(n_samples),
                "n_segregating": n_segregating,
            }

    return sfs_values


def calculate_fst(gt_array: allel.GenotypeArray, populations: Dict[str, List[int]]) -> pd.DataFrame:
    pop_names = list(populations.keys())
    fst_matrix = np.zeros((len(pop_names), len(pop_names)))

    logger.info("计算 FST: 群体数 %s", len(pop_names))

    for i, pop1 in enumerate(pop_names):
        for j, pop2 in enumerate(pop_names):
            if i < j:
                subpops = [populations[pop1], populations[pop2]]
                try:
                    a, b, c = allel.weir_cockerham_fst(gt_array, subpops)
                    denom = a + b + c
                    denom_sum = np.sum(denom)
                    if denom_sum == 0:
                        fst_value = 0.0
                    else:
                        fst_value = float(np.sum(a) / denom_sum)
                    fst_value = float(np.clip(fst_value, 0.0, 1.0))
                    fst_matrix[i, j] = fst_matrix[j, i] = fst_value
                except Exception as exc:
                    logger.warning("计算 %s vs %s 的 FST 时出错: %s", pop1, pop2, exc)
                    fst_matrix[i, j] = fst_matrix[j, i] = 0.0

    return pd.DataFrame(fst_matrix, index=pop_names, columns=pop_names)


def calculate_dxy(
    gt_array: allel.GenotypeArray,
    populations: Dict[str, List[int]],
    pos: np.ndarray,
) -> pd.DataFrame:
    pop_names = list(populations.keys())
    dxy_matrix = np.zeros((len(pop_names), len(pop_names)))

    logger.info("计算 Dxy: 群体数 %s", len(pop_names))

    for i, pop1 in enumerate(pop_names):
        for j, pop2 in enumerate(pop_names):
            if i < j:
                pop1_gt = gt_array.take(populations[pop1], axis=1)
                pop2_gt = gt_array.take(populations[pop2], axis=1)
                ac1 = pop1_gt.count_alleles()
                ac2 = pop2_gt.count_alleles()
                dxy = allel.sequence_divergence(pos, ac1, ac2)
                dxy_matrix[i, j] = dxy_matrix[j, i] = dxy

    return pd.DataFrame(dxy_matrix, index=pop_names, columns=pop_names)


def prepare_population_allele_counts(
    gt_array: allel.GenotypeArray,
    populations: Dict[str, List[int]],
) -> Dict[str, allel.AlleleCountsArray]:
    ac_by_pop = {}
    for pop_name, pop_indices in populations.items():
        pop_gt = gt_array.take(pop_indices, axis=1)
        ac_by_pop[pop_name] = pop_gt.count_alleles()
    return ac_by_pop


def bootstrap_diversity(
    ac_by_pop: Dict[str, allel.AlleleCountsArray],
    pos: np.ndarray,
    n_boot: int,
    seed: int | None = None,
) -> Dict[str, Dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    n_sites = int(pos.shape[0])
    if n_sites == 0:
        raise ValueError("无可用位点进行 bootstrap")

    boot_indices = rng.integers(0, n_sites, size=(n_boot, n_sites))
    results: Dict[str, Dict[str, np.ndarray]] = {}

    for pop in ac_by_pop:
        results[pop] = {
            "pi": np.zeros(n_boot, dtype=float),
            "theta_w": np.zeros(n_boot, dtype=float),
        }

    for b in range(n_boot):
        idx = boot_indices[b]
        pos_b = pos[idx]
        order = np.argsort(pos_b)
        pos_b = pos_b[order]
        for pop, ac in ac_by_pop.items():
            ac_b = ac[idx]
            ac_b = ac_b[order]
            results[pop]["pi"][b] = allel.sequence_diversity(pos_b, ac_b)
            results[pop]["theta_w"][b] = allel.watterson_theta(pos_b, ac_b)

    return results


def permutation_test_fst_dxy(
    gt_array: allel.GenotypeArray,
    populations: Dict[str, List[int]],
    pos: np.ndarray,
    n_perm: int,
    seed: int | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    pop_names = list(populations.keys())
    fst_p = np.zeros((len(pop_names), len(pop_names)))
    dxy_p = np.zeros((len(pop_names), len(pop_names)))

    for i, pop1 in enumerate(pop_names):
        for j, pop2 in enumerate(pop_names):
            if i < j:
                idx1 = populations[pop1]
                idx2 = populations[pop2]
                combined = np.array(idx1 + idx2)
                size1 = len(idx1)
                size2 = len(idx2)

                subpops = [idx1, idx2]
                a, b, c = allel.weir_cockerham_fst(gt_array, subpops)
                denom = a + b + c
                denom_sum = np.sum(denom)
                obs_fst = 0.0 if denom_sum == 0 else float(np.sum(a) / denom_sum)
                obs_fst = float(np.clip(obs_fst, 0.0, 1.0))

                pop1_gt = gt_array.take(idx1, axis=1)
                pop2_gt = gt_array.take(idx2, axis=1)
                ac1 = pop1_gt.count_alleles()
                ac2 = pop2_gt.count_alleles()
                obs_dxy = float(allel.sequence_divergence(pos, ac1, ac2))

                fst_null = np.zeros(n_perm, dtype=float)
                dxy_null = np.zeros(n_perm, dtype=float)

                for k in range(n_perm):
                    rng.shuffle(combined)
                    p1 = combined[:size1].tolist()
                    p2 = combined[size1 : size1 + size2].tolist()

                    a_k, b_k, c_k = allel.weir_cockerham_fst(gt_array, [p1, p2])
                    denom_k = a_k + b_k + c_k
                    denom_sum_k = np.sum(denom_k)
                    fst_k = 0.0 if denom_sum_k == 0 else float(np.sum(a_k) / denom_sum_k)
                    fst_null[k] = float(np.clip(fst_k, 0.0, 1.0))

                    pop1_k = gt_array.take(p1, axis=1)
                    pop2_k = gt_array.take(p2, axis=1)
                    ac1_k = pop1_k.count_alleles()
                    ac2_k = pop2_k.count_alleles()
                    dxy_null[k] = float(allel.sequence_divergence(pos, ac1_k, ac2_k))

                fst_p_val = (1.0 + np.sum(fst_null >= obs_fst)) / (n_perm + 1.0)
                dxy_p_val = (1.0 + np.sum(dxy_null >= obs_dxy)) / (n_perm + 1.0)

                fst_p[i, j] = fst_p[j, i] = fst_p_val
                dxy_p[i, j] = dxy_p[j, i] = dxy_p_val

    return (
        pd.DataFrame(fst_p, index=pop_names, columns=pop_names),
        pd.DataFrame(dxy_p, index=pop_names, columns=pop_names),
    )
