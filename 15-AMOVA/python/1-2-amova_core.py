"""
AMOVA（Analysis of Molecular Variance）核心算法库。

实现 Excoffier et al. (1992) Genetics 131:479-491 的等位基因频率方法，
适用于二等位基因单倍体 SNP 数据（线粒体 DNA）。

仅供 import 使用，无 CLI 入口。

算法概述：
  SSD 计算使用等位基因频率方法（避免构建 N×N 距离矩阵）：
    SSD_total       = Σ_l [ N * p_l * (1 - p_l) ]
    SSD_within      = Σ_l Σ_k [ n_k * p_kl * (1 - p_kl) ]
    SSD_among_pops  = SSD_within_group - SSD_within
    SSD_among_groups = SSD_total - SSD_within_group
"""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any

import numpy as np

log = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────
# SSD 计算
# ──────────────────────────────────────────────────────────────────────

def compute_ssd(
    pop_freq: np.ndarray,              # (n_pops, n_snps) float64/float32
    pop_sizes: np.ndarray,             # (n_pops,) int
    group_labels: np.ndarray | None,   # (n_pops,) 整数或字符串；None = 2级 AMOVA
) -> dict[str, float]:
    """
    计算 AMOVA 所需的各层 SSD（偏差平方和）。

    对于二等位基因单倍体数据：
      SSD_within = Σ_l Σ_k [ n_k * p_kl * (1 - p_kl) ]

    Args:
        pop_freq     — 各种群各位点 ALT 等位基因频率，(n_pops, n_snps)
        pop_sizes    — 各种群样本量
        group_labels — 各种群所属大组标签；None 时按 2级 AMOVA 处理

    Returns:
        含 ssd_total, ssd_within, ssd_among_pops, ssd_among_groups 的字典
    """
    pop_freq = np.asarray(pop_freq, dtype=np.float64)
    pop_sizes = np.asarray(pop_sizes, dtype=np.float64)
    n_total = float(pop_sizes.sum())

    # 全局等位基因频率（加权平均）
    p_global = (pop_freq * pop_sizes[:, None]).sum(axis=0) / n_total  # (n_snps,)

    # SSD_total = Σ_l [ N * p_l * (1 - p_l) ]
    ssd_total = float((n_total * p_global * (1.0 - p_global)).sum())

    # SSD_within = Σ_l Σ_k [ n_k * p_kl * (1 - p_kl) ]
    ssd_within = float((pop_sizes[:, None] * pop_freq * (1.0 - pop_freq)).sum())

    # 2级：所有种群视为一个大组
    if group_labels is None:
        ssd_within_group = ssd_total   # 整体与 total 相同（无大组层）
        ssd_among_groups = 0.0
        ssd_among_pops = ssd_total - ssd_within
        return {
            "ssd_total": ssd_total,
            "ssd_within": ssd_within,
            "ssd_among_pops": ssd_among_pops,
            "ssd_among_groups": ssd_among_groups,
            "ssd_within_group": ssd_within_group,
        }

    # 3级：计算各大组内等位基因频率并求 SSD_within_group
    group_arr = np.asarray(group_labels)
    unique_groups = np.unique(group_arr)
    ssd_within_group = 0.0
    for g in unique_groups:
        mask = group_arr == g
        n_g = float(pop_sizes[mask].sum())
        if n_g == 0:
            continue
        # 组内加权平均频率
        p_gl = (pop_freq[mask] * pop_sizes[mask, None]).sum(axis=0) / n_g
        ssd_within_group += float((n_g * p_gl * (1.0 - p_gl)).sum())

    ssd_among_pops = ssd_within_group - ssd_within
    ssd_among_groups = ssd_total - ssd_within_group

    return {
        "ssd_total": ssd_total,
        "ssd_within": ssd_within,
        "ssd_among_pops": ssd_among_pops,
        "ssd_among_groups": ssd_among_groups,
        "ssd_within_group": ssd_within_group,
    }


# ──────────────────────────────────────────────────────────────────────
# 方差分量系数（Excoffier 1992 Appendix）
# ──────────────────────────────────────────────────────────────────────

def compute_nc_coefficients(
    pop_sizes: np.ndarray,             # (n_pops,)
    group_labels: np.ndarray | None,   # (n_pops,) 或 None
) -> tuple[float, float]:
    """
    计算 AMOVA 方差分量系数 n̄_c 和 ñ_c（Excoffier 1992）。

    n̄_c：E[MSP] = Vc + n̄_c * Vb 中的系数
    ñ_c：E[MSG] = Vc + n̄_c * Vb + ñ_c * Va 中的系数

    2级 AMOVA（group_labels=None）仅使用 n̄_c；ñ_c 返回 np.nan。
    """
    pop_sizes = np.asarray(pop_sizes, dtype=np.float64)
    n_total = float(pop_sizes.sum())
    r = len(pop_sizes)      # 种群数

    if group_labels is None:
        # 2级：n̄_c = [N - Σ_k n²_k / N] / (r - 1)
        if r <= 1:
            return np.nan, np.nan
        nc_bar = (n_total - float((pop_sizes ** 2).sum()) / n_total) / (r - 1)
        return float(nc_bar), np.nan

    # 3级
    group_arr = np.asarray(group_labels)
    unique_groups = np.unique(group_arr)
    s = len(unique_groups)  # 大组数

    # Σ_γ (1/n_γ) * Σ_{α∈γ} n²_α
    inner_sum = 0.0
    for g in unique_groups:
        mask = group_arr == g
        n_g = float(pop_sizes[mask].sum())
        if n_g == 0:
            continue
        inner_sum += float((pop_sizes[mask] ** 2).sum()) / n_g

    # n̄_c = [N - inner_sum] / (r - s)
    df_p = r - s
    if df_p <= 0:
        nc_bar = np.nan
    else:
        nc_bar = (n_total - inner_sum) / df_p

    # ñ_c = (N - Σ_γ n²_γ / N) / (s - 1)
    # 该系数对应 E[MSG] 中 σ²_a 的系数，等价于组间等效样本量
    # 与 n̄_c 公式类比（在组层面进行同样的调整）
    n_g_vals = np.array([
        float(pop_sizes[group_arr == g].sum()) for g in unique_groups
    ])
    ng_sq_sum = float((n_g_vals ** 2).sum()) / n_total   # Σ_γ n²_γ / N

    df_g = s - 1
    if df_g <= 0:
        nc_tilde = np.nan
    else:
        nc_tilde = (n_total - ng_sq_sum) / df_g

    return float(nc_bar), float(nc_tilde)


# ──────────────────────────────────────────────────────────────────────
# 自由度
# ──────────────────────────────────────────────────────────────────────

def compute_degrees_of_freedom(
    pop_sizes: np.ndarray,
    group_labels: np.ndarray | None,
) -> dict[str, int]:
    """
    计算 AMOVA 各层自由度。

    Returns:
        df_among_groups, df_among_pops, df_within（2级时 df_among_groups=0）
    """
    pop_sizes = np.asarray(pop_sizes, dtype=np.int64)
    n_total = int(pop_sizes.sum())
    r = len(pop_sizes)

    if group_labels is None:
        return {
            "df_among_groups": 0,
            "df_among_pops": r - 1,
            "df_within": n_total - r,
        }

    group_arr = np.asarray(group_labels)
    s = len(np.unique(group_arr))
    return {
        "df_among_groups": s - 1,
        "df_among_pops": r - s,
        "df_within": n_total - r,
    }


# ──────────────────────────────────────────────────────────────────────
# 方差分量求解
# ──────────────────────────────────────────────────────────────────────

def solve_variance_components(
    ssd: dict[str, float],
    df: dict[str, int],
    nc_bar: float,
    nc_tilde: float,
    two_level: bool,
) -> tuple[float, float, float]:
    """
    从均方（MS）反解方差分量 Va、Vb、Vc。

    期望均方：
      E[MSW] = Vc
      E[MSP] = Vc + n̄_c * Vb
      E[MSG] = Vc + n̄_c * Vb + ñ_c * Va

    Args:
        ssd       — compute_ssd() 的返回值
        df        — compute_degrees_of_freedom() 的返回值
        nc_bar    — 系数 n̄_c
        nc_tilde  — 系数 ñ_c（2级时忽略）
        two_level — True=2级 AMOVA（无超级分组）

    Returns:
        (Va, Vb, Vc) — 3个方差分量（Va 为大组间，Vb 为组内种群间，Vc 为种群内）
    """
    df_g = df["df_among_groups"]
    df_p = df["df_among_pops"]
    df_w = df["df_within"]

    # 均方
    ms_w = ssd["ssd_within"] / df_w if df_w > 0 else 0.0
    ms_p = ssd["ssd_among_pops"] / df_p if df_p > 0 else 0.0

    vc = ms_w  # E[MSW] = Vc

    # Vb = (MSP - Vc) / n̄_c
    if np.isnan(nc_bar) or nc_bar == 0:
        vb = 0.0
    else:
        vb = (ms_p - vc) / nc_bar

    if two_level:
        return 0.0, vb, vc

    # Va = (MSG - MSP) / ñ_c
    ms_g = ssd["ssd_among_groups"] / df_g if df_g > 0 else 0.0
    if np.isnan(nc_tilde) or nc_tilde == 0:
        va = 0.0
    else:
        va = (ms_g - ms_p) / nc_tilde

    return va, vb, vc


# ──────────────────────────────────────────────────────────────────────
# F 统计量
# ──────────────────────────────────────────────────────────────────────

def compute_f_statistics(
    va: float,
    vb: float,
    vc: float,
    two_level: bool,
) -> dict[str, float]:
    """
    从方差分量计算 F 统计量（FST、FSC、FCT）。

    FST = (Va + Vb) / (Va + Vb + Vc)   — 总体种群分化
    FSC = Vb / (Vb + Vc)                — 大组内种群间分化
    FCT = Va / (Va + Vb + Vc)           — 大组间分化

    2级 AMOVA：FST = Vb / (Vb + Vc)，FSC/FCT 均为 NaN。
    """
    v_total = va + vb + vc
    if v_total == 0:
        return {"fst": np.nan, "fsc": np.nan, "fct": np.nan}

    fst = (va + vb) / v_total
    if two_level:
        return {"fst": float(fst), "fsc": np.nan, "fct": np.nan}

    fsc = vb / (vb + vc) if (vb + vc) != 0 else np.nan
    fct = va / v_total

    return {"fst": float(fst), "fsc": float(fsc), "fct": float(fct)}


# ──────────────────────────────────────────────────────────────────────
# 置换检验
# ──────────────────────────────────────────────────────────────────────

def _perm_chunk_worker(
    gt_T: np.ndarray,          # (N, n_snps) C-order，只读共享
    n_chunk: int,              # 本批次置换次数
    boundaries: np.ndarray,   # 种群边界
    pop_sizes: np.ndarray,    # (n_pops,) int64
    pop_sizes_f: np.ndarray,  # (n_pops,) float64
    ssd_total: float,
    df_within: int,
    df_among_pops: int,
    nc_bar: float,
    obs_fst: float,
    seed: int,                 # 本批次独立随机种子
) -> int:
    """单批次置换工作函数，供 ThreadPoolExecutor 并行调用。"""
    rng = np.random.default_rng(seed)
    n_total = int(pop_sizes.sum())
    n_pops = len(pop_sizes)
    n_snps = gt_T.shape[1]
    alt_by_pop = np.empty((n_pops, n_snps), dtype=np.int32)
    nv_by_pop  = np.empty((n_pops, n_snps), dtype=np.int32)
    count_ge = 0
    for _ in range(n_chunk):
        perm_idx = rng.permutation(n_total)
        for k in range(n_pops):
            b0, b1 = int(boundaries[k]), int(boundaries[k + 1])
            sub = gt_T[perm_idx[b0:b1], :]
            alt_by_pop[k] = (sub == 1).sum(axis=0, dtype=np.int32)
            nv_by_pop[k]  = (sub >= 0).sum(axis=0, dtype=np.int32)
        safe_nv = np.where(nv_by_pop > 0, nv_by_pop.astype(np.float64), 1.0)
        p_kl = alt_by_pop.astype(np.float64) / safe_nv
        ssd_w = float((pop_sizes_f[:, None] * p_kl * (1.0 - p_kl) * (nv_by_pop > 0)).sum())
        ssd_p = ssd_total - ssd_w
        ms_w = ssd_w / df_within if df_within > 0 else 0.0
        ms_p = ssd_p / df_among_pops if df_among_pops > 0 else 0.0
        vc = ms_w
        vb = (ms_p - vc) / nc_bar if not (np.isnan(nc_bar) or nc_bar == 0) else 0.0
        denom = vb + vc
        fst = vb / denom if denom > 0 else 0.0
        if fst >= obs_fst:
            count_ge += 1
    return count_ge


def permutation_test(
    gt_matrix: np.ndarray,             # (n_snps, N) int8，样本按种群顺序排列
    pop_sizes: np.ndarray,             # (n_pops,)
    group_labels: np.ndarray | None,
    obs_fst: float,
    n_perm: int,
    rng: np.random.Generator,
    two_level: bool,
    ssd_total: float,                  # 预计算的常数，避免重复计算
    df_within: int,
    df_among_pops: int,
    nc_bar: float,
    n_jobs: int = 1,                   # 并行线程数（1=串行）
) -> float:
    """
    个体水平置换检验：打乱个体-种群分配，保持种群大小不变，重计算 FST。

    P-value = (count_perm_FST ≥ obs_FST + 1) / (n_perm + 1)

    n_jobs > 1 时将 n_perm 均分给多个线程（ThreadPoolExecutor）并行执行。
    numpy C 层操作会释放 GIL，多线程可实现真并行。
    """
    pop_sizes = np.asarray(pop_sizes, dtype=np.int64)
    n_pops = len(pop_sizes)
    n_total = int(pop_sizes.sum())
    pop_sizes_f = pop_sizes.astype(np.float64)
    boundaries = np.concatenate([[0], np.cumsum(pop_sizes)]).astype(int)
    gt_T = np.ascontiguousarray(gt_matrix.T)   # (N, n_snps) C-order，行索引缓存友好

    # 生成各批次独立种子（在主线程内，确保可复现）
    n_workers = max(1, n_jobs)
    chunk_sizes = [n_perm // n_workers] * n_workers
    for i in range(n_perm % n_workers):
        chunk_sizes[i] += 1
    seeds = rng.integers(0, 2**31, size=n_workers)

    kwargs = dict(
        gt_T=gt_T,
        boundaries=boundaries,
        pop_sizes=pop_sizes,
        pop_sizes_f=pop_sizes_f,
        ssd_total=ssd_total,
        df_within=df_within,
        df_among_pops=df_among_pops,
        nc_bar=nc_bar,
        obs_fst=obs_fst,
    )

    if n_workers == 1:
        count_ge = _perm_chunk_worker(n_chunk=chunk_sizes[0], seed=int(seeds[0]), **kwargs)
    else:
        count_ge = 0
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = [
                pool.submit(_perm_chunk_worker,
                            n_chunk=chunk_sizes[i], seed=int(seeds[i]), **kwargs)
                for i in range(n_workers)
            ]
            for f in as_completed(futures):
                count_ge += f.result()

    p_value = (count_ge + 1) / (n_perm + 1)
    return float(p_value)


# ──────────────────────────────────────────────────────────────────────
# 顶层入口：完整 AMOVA 流程
# ──────────────────────────────────────────────────────────────────────

def run_amova(
    pop_freq: np.ndarray,              # (n_pops, n_snps)
    pop_sizes: np.ndarray,             # (n_pops,)
    group_labels: np.ndarray | None,   # (n_pops,) 或 None
    gt_matrix: np.ndarray,             # (n_snps, N) int8
    pop_names: list[str],              # 各种群名称（用于日志）
    n_perm: int,
    random_seed: int,
    n_jobs: int = 1,                   # 置换检验并行线程数
) -> dict[str, Any]:
    """
    执行完整 AMOVA 流程：SSD 计算 → 方差分量 → F 统计量 → 置换检验。

    Args:
        pop_freq    — (n_pops, n_snps) 等位基因频率矩阵
        pop_sizes   — (n_pops,) 各种群样本量
        group_labels — (n_pops,) 大组标签；None = 2级 AMOVA
        gt_matrix   — (n_snps, N) 基因型矩阵（-1=缺失）
        pop_names   — 种群名称列表
        n_perm      — 置换次数
        random_seed — 随机种子

    Returns:
        包含所有统计量的字典
    """
    two_level = group_labels is None
    pop_sizes = np.asarray(pop_sizes, dtype=np.int64)
    n_pops = len(pop_sizes)
    n_total = int(pop_sizes.sum())

    if group_labels is not None:
        group_arr = np.asarray(group_labels)
        n_groups = len(np.unique(group_arr))
    else:
        n_groups = 1

    log.info("运行 AMOVA：%d 种群，%d 大组，%d 样本，%s级结构",
             n_pops, n_groups, n_total, "2" if two_level else "3")

    # 1. SSD 计算
    ssd = compute_ssd(pop_freq, pop_sizes, group_labels)
    log.info("SSD — total=%.4f  within=%.4f  among_pops=%.4f  among_groups=%.4f",
             ssd["ssd_total"], ssd["ssd_within"],
             ssd["ssd_among_pops"], ssd["ssd_among_groups"])

    # 2. 自由度
    df = compute_degrees_of_freedom(pop_sizes, group_labels)
    log.info("自由度 — groups=%d  pops=%d  within=%d",
             df["df_among_groups"], df["df_among_pops"], df["df_within"])

    # 3. 方差分量系数
    nc_bar, nc_tilde = compute_nc_coefficients(pop_sizes, group_labels)
    log.info("系数 — n̄_c=%.4f  ñ_c=%s",
             nc_bar if not np.isnan(nc_bar) else float("nan"),
             f"{nc_tilde:.4f}" if not np.isnan(nc_tilde) else "N/A")

    # 4. 方差分量
    va, vb, vc = solve_variance_components(ssd, df, nc_bar, nc_tilde, two_level)
    v_total = va + vb + vc
    log.info("方差分量 — Va=%.6f  Vb=%.6f  Vc=%.6f  Total=%.6f",
             va, vb, vc, v_total)

    # 5. F 统计量
    f_stats = compute_f_statistics(va, vb, vc, two_level)
    log.info("F统计量 — FST=%.5f  FSC=%s  FCT=%s",
             f_stats["fst"],
             f"{f_stats['fsc']:.5f}" if not np.isnan(f_stats["fsc"]) else "N/A",
             f"{f_stats['fct']:.5f}" if not np.isnan(f_stats["fct"]) else "N/A")

    # 6. 百分比方差分量
    if v_total != 0:
        pct_among_groups = 100.0 * va / v_total
        pct_among_pops = 100.0 * vb / v_total
        pct_within = 100.0 * vc / v_total
    else:
        pct_among_groups = pct_among_pops = pct_within = np.nan

    # 7. 置换检验（基于 FST）
    log.info("开始置换检验（n=%d，线程数=%d）…", n_perm, n_jobs)
    rng = np.random.default_rng(random_seed)
    p_value = permutation_test(
        gt_matrix=gt_matrix,
        pop_sizes=pop_sizes,
        group_labels=group_labels,
        obs_fst=f_stats["fst"],
        n_perm=n_perm,
        rng=rng,
        two_level=two_level,
        ssd_total=ssd["ssd_total"],
        df_within=df["df_within"],
        df_among_pops=df["df_among_pops"],
        nc_bar=nc_bar,
        n_jobs=n_jobs,
    )
    log.info("置换检验 P-value = %.5f", p_value)

    return {
        # 基本信息
        "n_pops": n_pops,
        "n_groups": n_groups,
        "n_total": n_total,
        "two_level": two_level,
        # SSD
        "ssd_total": ssd["ssd_total"],
        "ssd_within": ssd["ssd_within"],
        "ssd_among_pops": ssd["ssd_among_pops"],
        "ssd_among_groups": ssd["ssd_among_groups"],
        # 自由度
        "df_among_groups": df["df_among_groups"],
        "df_among_pops": df["df_among_pops"],
        "df_within": df["df_within"],
        # 方差分量
        "va": va,
        "vb": vb,
        "vc": vc,
        # 百分比
        "pct_among_groups": pct_among_groups,
        "pct_among_pops": pct_among_pops,
        "pct_within": pct_within,
        # F 统计量
        "fst": f_stats["fst"],
        "fsc": f_stats["fsc"],
        "fct": f_stats["fct"],
        # 置换检验
        "p_value": p_value,
        "n_perm": n_perm,
    }
