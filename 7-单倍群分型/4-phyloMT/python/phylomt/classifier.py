"""单倍群分类算法。"""

from __future__ import annotations

import multiprocessing
from typing import TYPE_CHECKING

from .models import AnnotationEntry, ClassificationHit, SampleProfile, TreeBundle, TreeNode
from .mutations import base_token, mutation_position, position_covered, sum_weights

if TYPE_CHECKING:
    pass

# 进程池 initializer 使用的全局上下文（子进程中）
_worker_node_profiles: list[tuple[str, list[str], set[str], set[int], list[str]]] = []
_worker_node_order: dict[str, int] = {}
_worker_weights: dict[str, float] = {}
_worker_hotspots: set[str] = set()
_worker_annotation_db: dict[str, AnnotationEntry] = {}


def _init_worker(node_profiles, node_order, weights, hotspots, annotation_db):
    global _worker_node_profiles, _worker_node_order, _worker_weights
    global _worker_hotspots, _worker_annotation_db
    _worker_node_profiles = node_profiles
    _worker_node_order = node_order
    _worker_weights = weights
    _worker_hotspots = hotspots
    _worker_annotation_db = annotation_db


def _annotate_extra(token: str) -> str:
    """将单个 extra 变异分类为 hotspot / localPrivateMutation / globalPrivateMutation。"""
    if token in _worker_hotspots:
        return "hotspot"
    entry = _worker_annotation_db.get(token)
    if entry is None:
        return "globalPrivateMutation"
    if entry.in_phylotree:
        return "localPrivateMutation"
    return "globalPrivateMutation"


def _format_aac(token: str) -> str:
    """若变异有氨基酸变化注释，返回格式化字符串，否则返回空字符串。"""
    entry = _worker_annotation_db.get(token)
    if entry is None or not entry.aac:
        return ""
    return f"{token} [{entry.aac}| Codon {entry.codon_pos} | {entry.gene} ]"


def _confirm_back_mutations(
    path_back_muts: list[str],
    path_variant_positions: set[int],
    sample_positions: set[int],
    covered: "set[int] | None",
) -> list[str]:
    """检查路径上的返祖变异是否被样本确认（样本在该位置无变异）。"""
    confirmed: list[str] = []
    for bm_token in path_back_muts:
        bm_base = base_token(bm_token)
        pos = mutation_position(bm_base)
        if pos is None:
            continue
        if not position_covered(bm_base, covered):
            continue
        if pos in path_variant_positions:
            continue
        if pos not in sample_positions:
            confirmed.append(bm_base + "!")
    return confirmed


def _classify_one(args: tuple) -> tuple[str, list[ClassificationHit]]:
    """子进程执行单个样本的分类。"""

    sample_id, range_text, variants, covered_positions, hits = args
    sample_set = set(variants)
    sample_positions = {mutation_position(v) for v in variants if mutation_position(v) is not None}
    covered = covered_positions
    sample_weight = sum_weights(variants, _worker_weights)
    sample_hits: list[ClassificationHit] = []

    for name, path_variants, expected_set, path_variant_positions, path_back_muts in _worker_node_profiles:
        covered_expected = [v for v in path_variants if position_covered(v, covered)]
        covered_expected_set = set(covered_expected)
        found = sorted(sample_set & covered_expected_set, key=_variant_sort_key)
        missing = sorted(covered_expected_set - sample_set, key=_variant_sort_key)
        extra = sorted(sample_set - expected_set, key=_variant_sort_key)
        found_weight = sum_weights(found, _worker_weights)
        expected_weight = sum_weights(covered_expected, _worker_weights)
        quality = kulczynski(found_weight, sample_weight, expected_weight)

        confirmed_back = _confirm_back_mutations(path_back_muts, path_variant_positions, sample_positions, covered)

        annotated_extra = [(v, _annotate_extra(v)) for v in extra]
        aac_list = [s for v, _ in annotated_extra if (s := _format_aac(v))]

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
                confirmed_back_mutations=confirmed_back,
                annotated_extra_variants=annotated_extra,
                aac_in_remainings=aac_list,
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
    node_profiles = [
        (
            node.name,
            node.path_variants,
            set(node.path_variants),
            {mutation_position(v) for v in node.path_variants if mutation_position(v) is not None},
            node.path_back_mutations,
        )
        for node in nodes
    ]
    node_order = {node.name: index for index, node in enumerate(nodes)}
    weights = tree_bundle.weights
    hotspots = tree_bundle.hotspots
    annotation_db = tree_bundle.annotation_db

    tasks = [
        (s.sample_id, s.range_text, s.variants, s.covered_positions, hits)
        for s in samples
    ]

    if threads > 1 and len(samples) > 1:
        ctx = multiprocessing.get_context("fork")
        with ctx.Pool(
            processes=min(threads, len(samples)),
            initializer=_init_worker,
            initargs=(node_profiles, node_order, weights, hotspots, annotation_db),
        ) as pool:
            results_list = pool.map(_classify_one, tasks)
    else:
        _init_worker(node_profiles, node_order, weights, hotspots, annotation_db)
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
