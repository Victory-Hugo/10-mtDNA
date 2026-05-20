"""
读取步骤1输出的频率矩阵，遍历 scenarios.yaml 中的分析场景，
为每个场景独立运行 AMOVA 并将结果写入 TSV 文件。
"""

import argparse
import hashlib
import json
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


def _compute_scenario_key(
    step1_key: str,
    scenario: dict,
    permutation_n: int,
    scene_seed: int,
) -> str:
    """
    计算单个场景的缓存指纹。

    组成：Step1 指纹 + 场景定义（含分组文件路径/列名/subset 等全部字段）
          + 置换次数 + 本场景随机种子。
    任意一项变化均会导致缓存失效，触发重新计算。
    """
    src = json.dumps(
        {
            "step1_key": step1_key,
            "scenario": scenario,
            "permutation_n": permutation_n,
            "scene_seed": scene_seed,
        },
        sort_keys=True,
        ensure_ascii=False,
    )
    return hashlib.md5(src.encode()).hexdigest()


_SCENARIO_CACHE_FILENAME = ".amova_cache"


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
    step1_key: str = "",
) -> dict:
    """
    遍历所有场景，逐一运行 AMOVA 并保存结果。

    断点续跑：每个场景独立维护缓存键（.amova_cache），
    由 step1_key + 场景定义 + 置换参数共同决定。
    overwrite=False 时命中缓存则跳过；overwrite=True 时强制重算并刷新缓存。

    Args:
        freq_matrix_path — 步骤1产出的 .npz 文件
        pop_names_path   — 步骤1产出的种群名称 .txt 文件
        scenarios_path   — scenarios.yaml 路径
        base_dir         — 项目根目录（用于解析场景中的相对路径）
        output_dir       — 场景结果输出根目录
        permutation_n    — 置换检验次数
        random_seed      — 随机种子（各场景种子递增，保证独立性）
        overwrite        — 是否覆盖已有结果
        n_jobs           — 置换检验并行线程数
        step1_key        — Step1 输入指纹（由 pipeline.sh 传入）

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
        scene_seed = random_seed + idx
        log.info("═══ 场景 %d/%d: %s ═══", idx + 1, len(scenarios), label)

        out_dir = out_root / name
        result_path = out_dir / "amova_result.tsv"
        cache_file = out_dir / _SCENARIO_CACHE_FILENAME

        # ── 缓存命中检查 ──────────────────────────────────────────────────────
        scenario_key = _compute_scenario_key(step1_key, scenario, permutation_n, scene_seed)

        if (
            not overwrite
            and result_path.exists()
            and cache_file.exists()
            and cache_file.read_text().strip() == scenario_key
        ):
            log.info("缓存命中，跳过场景: %s", name)
            results_summary.append({"name": name, "status": "cache_hit"})
            continue

        if not overwrite and result_path.exists() and not (
            cache_file.exists() and cache_file.read_text().strip() == scenario_key
        ):
            log.info("缓存失效（场景定义或参数已变更），重新计算: %s", name)

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

        # 3. 运行 AMOVA
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

        # 4. 保存结果 TSV 并更新缓存键
        df_result = pd.DataFrame([result])
        df_result.to_csv(result_path, sep="\t", index=False)
        cache_file.write_text(scenario_key)
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
    p.add_argument("--step1-key", default="",
                   help="Step1 输入指纹（由 pipeline.sh 传入，用于场景缓存）")
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
        step1_key=args.step1_key,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
