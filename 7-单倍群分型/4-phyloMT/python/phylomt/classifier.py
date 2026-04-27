"""单倍群分类算法。"""

from __future__ import annotations

import multiprocessing
from typing import TYPE_CHECKING

from .models import ClassificationHit, SampleProfile, TreeBundle, TreeNode
from .mutations import position_covered, sum_weights

if TYPE_CHECKING:
    pass

# 进程池 initializer 使用的全局上下文（子进程中）
_worker_node_profiles: list[tuple[str, list[str]]] = []
_worker_node_order: dict[str, int] = {}
_worker_weights: dict[str, float] = {}


def _init_worker(node_profiles, node_order, weights):
    global _worker_node_profiles, _worker_node_order, _worker_weights
    _worker_node_profiles = node_profiles
    _worker_node_order = node_order
    _worker_weights = weights


def _classify_one(args: tuple) -> tuple[str, list[ClassificationHit]]:
    """子进程执行单个样本的分类。"""

    sample_id, range_text, variants, covered_positions, hits = args
    sample_set = set(variants)
    covered = covered_positions
    sample_weight = sum_weights(variants, _worker_weights)
    sample_hits: list[ClassificationHit] = []

    for name, path_variants in _worker_node_profiles:
        covered_expected = [v for v in path_variants if position_covered(v, covered)]
        covered_expected_set = set(covered_expected)
        expected_set = set(path_variants)
        found = sorted(sample_set & covered_expected_set, key=_variant_sort_key)
        missing = sorted(covered_expected_set - sample_set, key=_variant_sort_key)
        extra = sorted(sample_set - expected_set, key=_variant_sort_key)
        found_weight = sum_weights(found, _worker_weights)
        expected_weight = sum_weights(covered_expected, _worker_weights)
        quality = kulczynski(found_weight, sample_weight, expected_weight)
        sample_hits.append(
            ClassificationHit(
                sample_id=sample_id,
                haplogroup=name,
                rank=0,
                quality=quality,
                range_text=range_text,
                found_variants=found,
                missing_variants=missing,
                extra_variants=extra,
                found_weight=found_weight,
                expected_weight=expected_weight,
                sample_weight=sample_weight,
            )
        )

    ranked = sorted(sample_hits, key=lambda h: (-h.quality, _worker_node_order[h.haplogroup]))[:hits]
    for index, hit in enumerate(ranked, start=1):
        hit.rank = index
    return sample_id, ranked


def classify_samples(
    samples: list[SampleProfile],
    tree_bundle: TreeBundle,
    hits: int,
    threads: int = 1,
) -> dict[str, list[ClassificationHit]]:
    """对多个样本执行分类，支持多进程并行。"""

    nodes = flatten_nodes(tree_bundle.root)
    node_profiles = [(node.name, node.path_variants) for node in nodes]
    node_order = {node.name: index for index, node in enumerate(nodes)}
    weights = tree_bundle.weights

    tasks = [
        (s.sample_id, s.range_text, s.variants, s.covered_positions, hits)
        for s in samples
    ]

    if threads > 1 and len(samples) > 1:
        ctx = multiprocessing.get_context("fork")
        with ctx.Pool(
            processes=min(threads, len(samples)),
            initializer=_init_worker,
            initargs=(node_profiles, node_order, weights),
        ) as pool:
            results_list = pool.map(_classify_one, tasks)
    else:
        _init_worker(node_profiles, node_order, weights)
        results_list = [_classify_one(t) for t in tasks]

    return {sid: hits_list for sid, hits_list in results_list}


def flatten_nodes(root: TreeNode) -> list[TreeNode]:
    """展平树节点。"""

    nodes = [root]
    for child in root.children:
        nodes.extend(flatten_nodes(child))
    return nodes


def kulczynski(found_weight: float, sample_weight: float, expected_weight: float) -> float:
    """Kulczynski 评分。"""

    sample_score = found_weight / sample_weight if sample_weight > 0 else 0.0
    expected_score = found_weight / expected_weight if expected_weight > 0 else 1.0
    return round(0.5 * sample_score + 0.5 * expected_score, 10)


def _variant_sort_key(token: str) -> tuple[int, str]:
    digits = ""
    for character in token:
        if character.isdigit():
            digits += character
        else:
            break
    return (int(digits) if digits else 0, token)
