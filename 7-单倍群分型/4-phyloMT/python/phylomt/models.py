"""核心数据结构。"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class SampleProfile:
    """标准化后的样本变异档案。"""

    sample_id: str
    range_text: str
    variants: list[str]
    source_path: Path
    ignored_positions: set[int] = field(default_factory=set)
    covered_positions: "set[int] | None" = field(default=None)


@dataclass(slots=True)
class AnnotationEntry:
    """单个变异的注释信息。"""

    category: str
    in_phylotree: bool
    aac: str
    codon_pos: str
    gene: str


@dataclass(slots=True)
class TreeNode:
    """树节点。"""

    name: str
    own_variants: list[str]
    children: list["TreeNode"] = field(default_factory=list)
    path_variants: list[str] = field(default_factory=list)
    path_back_mutations: list[str] = field(default_factory=list)


@dataclass(slots=True)
class TreeBundle:
    """树资源包。"""

    tree_id: str
    version: str
    tree_dir: Path
    reference_name: str
    reference_sequence: str
    root: TreeNode
    weights: dict[str, float]
    rules: list[tuple[list[str], list[str]]]
    hotspots: set[str] = field(default_factory=set)
    annotation_db: "dict[str, AnnotationEntry]" = field(default_factory=dict)

    @property
    def full_id(self) -> str:
        return f"{self.tree_id}@{self.version}"


@dataclass(slots=True)
class ClassificationHit:
    """单个候选单倍群命中。"""

    sample_id: str
    haplogroup: str
    rank: int
    quality: float
    range_text: str
    found_variants: list[str]
    missing_variants: list[str]
    extra_variants: list[str]
    found_weight: float
    expected_weight: float
    sample_weight: float
    confirmed_back_mutations: list[str] = field(default_factory=list)
    annotated_extra_variants: "list[tuple[str, str]]" = field(default_factory=list)
    aac_in_remainings: list[str] = field(default_factory=list)
