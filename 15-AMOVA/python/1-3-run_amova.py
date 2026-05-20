"""
读取步骤1输出的频率矩阵，遍历 scenarios.yaml 中的分析场景，
为每个场景独立运行 AMOVA 并将结果写入 TSV 文件。
"""

import argparse
import logging
from pathlib import Path

import importlib.util

import numpy as np
import pandas as pd
import yaml

# 本项目算法库：文件名含数字前缀，需用 importlib 加载
_core_path = Path(__file__).parent / "1-2-amova_core.py"
_spec = importlib.util.spec_from_file_location("amova_core", _core_path)
_core = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_core)

log = logging.getLogger(__name__)

_AMOVA_CORE = _core   # 允许测试时注入 mock


# ──────────────────────────────────────────────────────────────────────
# 辅助函数
# ──────────────────────────────────────────────────────────────────────

def _load_npz(freq_matrix_path: str) -> dict[str, np.ndarray]:
    """加载步骤1产出的 .npz 文件，返回各数组的字典。"""
    data = np.load(freq_matrix_path)
    return {
        "freq_matrix": data["freq_matrix"].astype(np.float64),
        "gt_matrix": data["gt_matrix"],
        "pop_sizes": data["pop_sizes"],
        "snp_positions": data["snp_positions"],
    }


def _load_pop_names(pop_names_path: str) -> list[str]:
    """加载步骤1产出的种群名称文件。"""
    with open(pop_names_path, encoding="utf-8") as fh:
        names = [line.strip() for line in fh if line.strip()]
    return names


def _load_scenarios(scenarios_path: str) -> list[dict]:
    """解析 scenarios.yaml，返回场景列表。"""
    with open(scenarios_path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    scenarios = cfg.get("scenarios", [])
    if not scenarios:
        raise ValueError(f"scenarios.yaml 中未找到任何场景: {scenarios_path}")
    return scenarios


def _build_group_labels(
    pop_names: list[str],
    scenario: dict,
    base_dir: str,
) -> np.ndarray | None:
    """
    根据场景配置构建大组标签数组（aligned with pop_names 顺序）。

    若 group_col 为 None，返回 None（2级 AMOVA）。
    若某个种群不在场景的分组文件中，该种群将被标记为 None 并在后续过滤。
    """
    group_col = scenario.get("group_col")
    if not group_col:
        return None

    group_file = str(Path(base_dir) / scenario["group_file"])
    id_col = scenario.get("id_col", "ID")
    pop_col = scenario.get("pop_col", "Group_small")
    sep = "\t" if Path(group_file).suffix in (".tsv", ".txt") else ","

    df = pd.read_csv(group_file, sep=sep, dtype=str)
    if pop_col not in df.columns or group_col not in df.columns:
        raise ValueError(
            f"分组文件 {group_file} 缺少列: {pop_col} 或 {group_col}"
        )

    # 每个种群只需要对应的大组标签（一行即可）
    pop_to_group = (
        df[[pop_col, group_col]]
        .dropna()
        .drop_duplicates(subset=pop_col)
        .set_index(pop_col)[group_col]
        .to_dict()
    )

    labels = []
    for pname in pop_names:
        labels.append(pop_to_group.get(pname, None))

    return np.array(labels, dtype=object)


def _apply_subset_filter(
    scenario: dict,
    pop_names: list[str],
    pop_sizes: np.ndarray,
    freq_matrix: np.ndarray,
    gt_matrix: np.ndarray,
    group_labels: np.ndarray | None,
    base_dir: str,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    """
    按场景 subset_filter 配置过滤种群。

    subset_col / subset_vals 指定保留哪些大组的种群（如只保留 Manchu）。
    同时去除 group_labels 中为 None（无大组映射）的种群。
    """
    group_col = scenario.get("group_col")
    subset_col = scenario.get("subset_col")
    subset_vals = scenario.get("subset_vals")

    # 构建过滤掩码（默认全部保留）
    keep_mask = np.ones(len(pop_names), dtype=bool)

    # 去除无大组映射的种群
    if group_labels is not None:
        for i, lbl in enumerate(group_labels):
            if lbl is None:
                keep_mask[i] = False

    # 子集过滤
    if subset_col and subset_vals and group_labels is not None:
        group_file = str(Path(base_dir) / scenario["group_file"])
        pop_col = scenario.get("pop_col", "Group_small")
        sep = "\t" if Path(group_file).suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(group_file, sep=sep, dtype=str)
        pop_to_subset = (
            df[[pop_col, subset_col]]
            .dropna()
            .drop_duplicates(subset=pop_col)
            .set_index(pop_col)[subset_col]
            .to_dict()
        )
        subset_set = set(str(v) for v in subset_vals)
        for i, pname in enumerate(pop_names):
            if str(pop_to_subset.get(pname, "")) not in subset_set:
                keep_mask[i] = False

    # 应用掩码
    keep_idx = np.where(keep_mask)[0]
    if len(keep_idx) == 0:
        raise ValueError("过滤后无剩余种群，请检查 subset_filter 配置")

    # 过滤 pop_names、pop_sizes、freq_matrix
    new_pop_names = [pop_names[i] for i in keep_idx]
    new_pop_sizes = pop_sizes[keep_idx]
    new_freq_matrix = freq_matrix[keep_idx, :]
    new_group_labels = group_labels[keep_idx] if group_labels is not None else None

    # 过滤 gt_matrix（列按种群顺序排列，需要根据 pop_sizes 重建样本边界）
    old_boundaries = np.concatenate([[0], np.cumsum(pop_sizes)]).astype(int)
    col_indices = np.concatenate([
        np.arange(old_boundaries[i], old_boundaries[i + 1])
        for i in keep_idx
    ])
    new_gt_matrix = gt_matrix[:, col_indices]

    return new_pop_names, new_pop_sizes, new_freq_matrix, new_gt_matrix, new_group_labels


# ──────────────────────────────────────────────────────────────────────
# 核心业务函数
# ──────────────────────────────────────────────────────────────────────

def run(
    freq_matrix_path: str,
    pop_names_path: str,
    scenarios_path: str,
    base_dir: str,
    output_dir: str,
    permutation_n: int,
    random_seed: int,
    overwrite: bool = True,
    n_jobs: int = 1,
) -> dict:
    """
    遍历所有场景，逐一运行 AMOVA 并保存结果。

    Args:
        freq_matrix_path — 步骤1产出的 .npz 文件
        pop_names_path   — 步骤1产出的种群名称 .txt 文件
        scenarios_path   — scenarios.yaml 路径
        base_dir         — 项目根目录（用于解析场景中的相对路径）
        output_dir       — 场景结果输出根目录
        permutation_n    — 置换检验次数
        random_seed      — 随机种子（各场景种子递增，保证独立性）
        overwrite        — 是否覆盖已有结果

    Returns:
        {"n_scenarios": int, "results": [per-scenario-dict]}
    """
    log.info("加载频率矩阵: %s", freq_matrix_path)
    npz = _load_npz(freq_matrix_path)
    base_pop_names = _load_pop_names(pop_names_path)

    if len(base_pop_names) != npz["freq_matrix"].shape[0]:
        raise ValueError(
            f"种群数量不一致: pop_names.txt 有 {len(base_pop_names)} 行，"
            f"freq_matrix 有 {npz['freq_matrix'].shape[0]} 列"
        )

    log.info("加载分析场景: %s", scenarios_path)
    scenarios = _load_scenarios(scenarios_path)
    log.info("共 %d 个场景", len(scenarios))

    out_root = Path(output_dir)
    results_summary = []

    for idx, scenario in enumerate(scenarios):
        name = scenario["name"]
        label = scenario.get("label", name)
        log.info("═══ 场景 %d/%d: %s ═══", idx + 1, len(scenarios), label)

        out_dir = out_root / name
        result_path = out_dir / "amova_result.tsv"

        if result_path.exists() and not overwrite:
            log.info("跳过（已存在）: %s", result_path)
            results_summary.append({"name": name, "status": "skipped"})
            continue

        out_dir.mkdir(parents=True, exist_ok=True)

        # 1. 构建大组标签
        group_labels = _build_group_labels(base_pop_names, scenario, base_dir)

        # 2. 子集过滤
        pop_names, pop_sizes, freq_matrix, gt_matrix, group_labels = _apply_subset_filter(
            scenario=scenario,
            pop_names=base_pop_names,
            pop_sizes=npz["pop_sizes"].copy(),
            freq_matrix=npz["freq_matrix"].copy(),
            gt_matrix=npz["gt_matrix"],
            group_labels=group_labels,
            base_dir=base_dir,
        )
        log.info("过滤后种群数: %d，样本数: %d", len(pop_names), int(pop_sizes.sum()))

        # 3. 运行 AMOVA（每个场景使用不同的种子保证独立性）
        scene_seed = random_seed + idx
        result = _AMOVA_CORE.run_amova(
            pop_freq=freq_matrix,
            pop_sizes=pop_sizes,
            group_labels=group_labels,
            gt_matrix=gt_matrix,
            pop_names=pop_names,
            n_perm=permutation_n,
            random_seed=scene_seed,
            n_jobs=n_jobs,
        )
        result["name"] = name
        result["label"] = label

        # 4. 保存结果 TSV
        df_result = pd.DataFrame([result])
        df_result.to_csv(result_path, sep="\t", index=False)
        log.info("已保存结果: %s", result_path)

        results_summary.append({"name": name, "status": "completed", "fst": result["fst"]})

    log.info("所有场景处理完成")
    return {"n_scenarios": len(scenarios), "results": results_summary}


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="遍历场景配置，为每个场景运行 AMOVA 并保存结果"
    )
    p.add_argument("--freq-matrix", required=True,
                   help="步骤1产出的 .npz 频率矩阵文件路径")
    p.add_argument("--pop-names", required=True,
                   help="步骤1产出的种群名称 .txt 文件路径")
    p.add_argument("--scenarios", required=True,
                   help="scenarios.yaml 文件路径")
    p.add_argument("--base-dir", required=True,
                   help="项目根目录（用于解析场景中的相对路径）")
    p.add_argument("--output-dir", required=True,
                   help="场景结果输出目录（每个场景一个子目录）")
    p.add_argument("--permutation-n", type=int, default=1000,
                   help="置换检验次数（默认: 1000）")
    p.add_argument("--random-seed", type=int, default=42,
                   help="随机种子（默认: 42）")
    p.add_argument("--overwrite", type=int, default=1,
                   help="是否覆盖已有结果（1=是，0=否，默认: 1）")
    p.add_argument("--n-jobs", type=int, default=1,
                   help="置换检验并行线程数（默认: 1，建议设为 CPU 核心数）")
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    args = build_parser().parse_args(argv)
    run(
        freq_matrix_path=args.freq_matrix,
        pop_names_path=args.pop_names,
        scenarios_path=args.scenarios,
        base_dir=args.base_dir,
        output_dir=args.output_dir,
        permutation_n=args.permutation_n,
        random_seed=args.random_seed,
        overwrite=bool(args.overwrite),
        n_jobs=args.n_jobs,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
