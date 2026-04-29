"""树资源扫描与树解析。"""

from __future__ import annotations

import gzip
import json
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from .models import AnnotationEntry, TreeBundle, TreeNode
from .mutations import apply_variants_to_profile, clean_token, is_back_mutation


PROJECT_DIR = Path(__file__).resolve().parents[2]
DEFAULT_TREES_DIR = PROJECT_DIR / "data" / "trees"
REQUIRED_TREE_METADATA_FIELDS = ("id", "version", "tree", "weights", "fasta", "alignmentRules")


def load_tree_metadata(tree_dir: Path) -> dict[str, str]:
    """读取树目录中的元数据。"""

    metadata: dict[str, str] = {}
    metadata_path = tree_dir / "tree.yaml"
    with metadata_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if line[:1].isspace() or stripped.startswith("- "):
                continue
            key, separator, value = stripped.partition(":")
            if not separator:
                continue
            value = value.strip()
            if not value:
                continue
            if value[:1] == value[-1:] and value[:1] in {'"', "'"}:
                value = value[1:-1]
            metadata[key] = value

    missing = [field for field in REQUIRED_TREE_METADATA_FIELDS if field not in metadata]
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(f"{metadata_path} 缺少必要字段: {missing_text}")
    return metadata


def list_installed_trees(trees_dir: Path) -> list[dict[str, str]]:
    """扫描目录中的已安装树。"""

    if not trees_dir.exists():
        return []

    entries: list[dict[str, str]] = []
    for tree_root in sorted(path for path in trees_dir.iterdir() if path.is_dir()):
        for version_dir in sorted(path for path in tree_root.iterdir() if path.is_dir()):
            if not (version_dir / "tree.yaml").exists():
                continue
            metadata = load_tree_metadata(version_dir)
            tree_id = tree_root.name
            version = version_dir.name
            entries.append(
                {
                    "id": tree_id,
                    "version": version,
                    "full_id": f"{tree_id}@{version}",
                    "name": str(metadata.get("name", "")),
                    "category": str(metadata.get("category", "")),
                    "path": str(version_dir),
                }
            )
    return entries


def load_tree_bundle(tree_dir: Path) -> TreeBundle:
    """从树目录加载树资源。"""

    metadata = load_tree_metadata(tree_dir)

    reference_record = next(SeqIO.parse(tree_dir / metadata["fasta"], "fasta"))
    reference_sequence = str(reference_record.seq).upper()
    tree_file = tree_dir / metadata["tree"]
    if tree_file.suffix == ".json":
        root = parse_tree_json(tree_file)
        fill_path_variants(root, [])
        # 若 tree.yaml 指定了 rsrsFasta，则将 RSRS 坐标系的路径变异转换为参考序列坐标系
        rsrs_fasta_name = metadata.get("rsrsFasta", "")
        if rsrs_fasta_name:
            rsrs_record = next(SeqIO.parse(tree_dir / rsrs_fasta_name, "fasta"))
            rsrs_sequence = str(rsrs_record.seq).upper()
            _convert_tree_path_variants_to_ref(root, rsrs_sequence, reference_sequence)
    else:
        root = parse_tree_xml(tree_file)
        fill_path_variants(root, [])
    weights = parse_weights(tree_dir / metadata["weights"])
    rules = parse_rules(tree_dir / metadata["alignmentRules"])
    hotspots = load_hotspots(tree_dir)
    annotation_db = load_annotation_db(tree_dir)

    return TreeBundle(
        tree_id=str(tree_dir.parent.name),
        version=str(tree_dir.name),
        tree_dir=tree_dir,
        reference_name=reference_record.id,
        reference_sequence=reference_sequence,
        root=root,
        weights=weights,
        rules=rules,
        hotspots=hotspots,
        annotation_db=annotation_db,
    )


_JSON_RANGE_DEL_RE = re.compile(r"^(\d+)-(\d+)d$", re.IGNORECASE)
_JSON_REF_PREFIX_RE = re.compile(r"^[A-Za-z](\d+.*)$")


def _normalize_json_variants(raw_token: str) -> list[str]:
    """将新版 JSON 变异格式规范化为旧版 XML 兼容格式（位置+碱基）。

    规则：
    - (G207A)   括号可选变异 → 跳过
    - reserved  占位符       → 跳过
    - 105-110d  范围删除     → 展开为 105d 106d … 110d
    - C498d     有参考碱基的删除 → 498d
    - A235G     有参考碱基替换  → 235G
    - T16311C!  有参考碱基返祖  → 16311C!
    - C152T!!   双返祖         → 152T!!
    - A16227c   小写碱基       → 16227C（大写）
    - 15883A    已是旧格式     → 保持原样
    - 5899.XC   插入           → 保持原样
    """
    token = raw_token.strip()
    if not token or token == "reserved":
        return []

    # 括号可选/复发性变异 → 跳过
    if token.startswith("(") and token.endswith(")"):
        return []

    # 范围删除 105-110d → [105d, 106d, …, 110d]
    m = _JSON_RANGE_DEL_RE.match(token)
    if m:
        start, end = int(m.group(1)), int(m.group(2))
        return [f"{pos}d" for pos in range(start, end + 1)]

    # 有参考碱基前缀的变异（单个字母开头）
    if token[0].isalpha():
        rest = token[1:]  # 去掉参考碱基

        # 双返祖 C152T!! → 152T!!
        if rest.endswith("!!"):
            body = rest[:-2]
            if body and body[-1].lower() in "acgt":
                body = body[:-1] + body[-1].upper()
            return [body + "!!"]

        # 单返祖 T16311C! → 16311C!
        if rest.endswith("!"):
            body = rest[:-1]
            if body and body[-1].lower() in "acgt":
                body = body[:-1] + body[-1].upper()
            return [body + "!"]

        # 有参考碱基的删除 C498d → 498d（尾部 d 不是碱基，保留小写）
        if rest and rest[-1].lower() not in "acgt":
            return [rest]

        # 普通替换或小写碱基 A235G → 235G, A16227c → 16227C
        if rest and rest[-1].isalpha():
            return [rest[:-1] + rest[-1].upper()]
        return [rest]

    # 已是旧格式（数字开头）：直接交给 clean_token 标准化
    return [clean_token(token)] if clean_token(token) else []


def parse_tree_json(json_path: Path) -> TreeNode:
    """解析新版 JSON 格式系统发育树。"""

    with json_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    haplogroups: dict = data["haplogroups"]

    # 构建 parent → children 映射
    children_map: dict[str, list[str]] = defaultdict(list)
    for name, info in haplogroups.items():
        parent = info.get("parent", "")
        children_map[parent].append(name)

    # 找根节点（parent 为空字符串）
    root_candidates = [
        name for name, info in haplogroups.items() if not info.get("parent")
    ]
    if not root_candidates:
        raise ValueError(f"{json_path} 中无根节点（parent 为空）。")
    root_name = root_candidates[0]

    def build_node(name: str) -> TreeNode:
        info = haplogroups[name]
        raw_mutations = info.get("mutations", "")
        own_variants: list[str] = []
        for raw_tok in raw_mutations.split():
            own_variants.extend(_normalize_json_variants(raw_tok))
        child_nodes = [build_node(c) for c in sorted(children_map.get(name, []))]
        return TreeNode(name=name, own_variants=own_variants, children=child_nodes)

    return build_node(root_name)


def _convert_tree_path_variants_to_ref(
    root: TreeNode, rsrs_seq: str, ref_seq: str
) -> None:
    """将 JSON 树所有节点的 path_variants 从 RSRS 坐标转换为目标参考序列坐标。

    前提：fill_path_variants 已运行，path_variants 中不含返祖标记（! 已被消解）。
    转换规则：
    - 对于 RSRS 与 ref_seq 一致的位点：path_variants 中的该变异不变。
    - 对于 RSRS 与 ref_seq 不一致的位点（共 N 个）：
        * 从 RSRS path_variants 中提取该节点在该位点的碱基（若无则取 RSRS 祖先碱基）
        * 若该碱基与 ref_seq（rCRS）不同 → 加入 ref 坐标变异
        * 若相同 → 不加（已是参考碱基）
    - 插入（含 .）、删除（以 d 结尾）：
        * 若位点不在差异集合中 → 保留
        * 若位点在差异集合中 → 按上述规则处理（通常很少见）
    """
    from .mutations import mutation_position

    # 构建 RSRS-ref 差异位点字典：{1-based pos: (ref_base, rsrs_base)}
    diff_map: dict[int, tuple[str, str]] = {}
    for i, (rb, sb) in enumerate(zip(ref_seq, rsrs_seq)):
        if rb != sb and rb != "N" and sb != "N":
            diff_map[i + 1] = (rb, sb)

    def convert_node(node: TreeNode) -> None:
        # 建立 RSRS 路径变异位点→碱基映射（仅 SNP 类型）
        pos_allele: dict[int, str] = {}
        non_snp: list[str] = []
        for tok in node.path_variants:
            pos = mutation_position(tok)
            if pos is None or "." in tok or tok.rstrip("!").endswith("d"):
                non_snp.append(tok)
            else:
                pos_allele[pos] = tok.lstrip("0123456789")  # 提取碱基部分

        rcrs_path: list[str] = []

        # 处理 RSRS-ref 差异位点
        for pos, (ref_base, rsrs_base) in diff_map.items():
            node_base = pos_allele.get(pos, rsrs_base)
            if node_base != ref_base:
                rcrs_path.append(f"{pos}{node_base}")

        # 保留非差异位点的 SNP 变异
        for tok in node.path_variants:
            pos = mutation_position(tok)
            if pos is not None and pos not in diff_map and "." not in tok and not tok.rstrip("!").endswith("d"):
                rcrs_path.append(tok)

        # 保留插入/删除变异（直接保留，位点与参考系无关）
        rcrs_path.extend(non_snp)

        node.path_variants = rcrs_path

        for child in node.children:
            convert_node(child)

    convert_node(root)


def parse_tree_xml(tree_path: Path) -> TreeNode:
    """解析 HaploGrep 树 XML。"""

    document_root = ET.parse(tree_path).getroot()
    first_node = document_root.find("haplogroup")
    if first_node is None:
        raise ValueError(f"{tree_path} 中没有 haplogroup 节点。")
    return parse_tree_node(first_node)


def parse_tree_node(element: ET.Element) -> TreeNode:
    """递归解析树节点。"""

    details = element.find("details")
    own_variants = [clean_token(poly.text or "") for poly in details.findall("poly")] if details is not None else []
    children = [parse_tree_node(child) for child in element.findall("haplogroup")]
    return TreeNode(
        name=element.attrib["name"],
        own_variants=[token for token in own_variants if token],
        children=children,
    )


def fill_path_variants(
    node: TreeNode,
    parent_variants: list[str],
    parent_back_muts: "list[str] | None" = None,
) -> None:
    """为每个节点预计算路径累计变异，同时追踪路径上的返祖变异。"""

    if parent_back_muts is None:
        parent_back_muts = []
    node.path_variants = apply_variants_to_profile(parent_variants, node.own_variants)
    own_back = [clean_token(t) for t in node.own_variants if is_back_mutation(clean_token(t))]
    node.path_back_mutations = parent_back_muts + own_back
    for child in node.children:
        fill_path_variants(child, node.path_variants, node.path_back_mutations)


def load_hotspots(tree_dir: Path) -> set[str]:
    """从 tree.yaml 读取 hotspot 变异列表。"""

    hotspots: set[str] = set()
    metadata_path = tree_dir / "tree.yaml"
    if not metadata_path.exists():
        return hotspots

    in_hotspots = False
    with metadata_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if stripped == "hotspots:":
                in_hotspots = True
                continue
            if in_hotspots:
                if stripped.startswith("- "):
                    item = stripped[2:].strip().strip("\"'")
                    if item:
                        hotspots.add(item)
                elif stripped and not line[:1].isspace():
                    in_hotspots = False
    return hotspots


def load_annotation_db(tree_dir: Path) -> "dict[str, AnnotationEntry]":
    """从 annotations/rCRS_annotation_*.txt.gz 加载变异注释数据库。"""

    annotation_dir = tree_dir / "annotations"
    if not annotation_dir.exists():
        return {}

    annotation_files = sorted(annotation_dir.glob("rCRS_annotation_*.txt.gz"))
    if not annotation_files:
        return {}

    db: dict[str, AnnotationEntry] = {}
    with gzip.open(annotation_files[0], "rt", encoding="utf-8") as handle:
        header_line = handle.readline().rstrip("\n")
        headers = header_line.split("\t")
        col = {name: i for i, name in enumerate(headers)}

        mut_col = col.get("Mutation", 0)
        cat_col = col.get("Category", 6)
        phylo_col = col.get("Phylotree17_haplogroups", 7)
        aac_col = col.get("AAC", 12)
        codon_col = col.get("CodonPosition", 13)
        locus_col = col.get("Maplocus", 5)

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= mut_col:
                continue
            mutation = fields[mut_col]
            if not mutation:
                continue

            def _get(idx: int) -> str:
                return fields[idx] if idx < len(fields) else ""

            category = _get(cat_col)
            phylo = _get(phylo_col)
            aac = _get(aac_col)
            codon_pos = _get(codon_col)
            maplocus = _get(locus_col)
            gene = maplocus[3:] if maplocus.startswith("MT-") else maplocus

            db[mutation] = AnnotationEntry(
                category=category,
                in_phylotree=bool(phylo.strip()),
                aac=aac,
                codon_pos=codon_pos,
                gene=gene,
            )
    return db


def parse_weights(weights_path: Path) -> dict[str, float]:
    """读取权重文件。"""

    weights: dict[str, float] = {}
    with weights_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            token = clean_token(fields[0])
            try:
                weights[token] = float(fields[1])
            except ValueError:
                continue
    return weights


def parse_rules(rules_path: Path) -> list[tuple[list[str], list[str]]]:
    """读取 FASTA 对齐规则。"""

    rules: list[tuple[list[str], list[str]]] = []
    with rules_path.open("r", encoding="utf-8") as handle:
        header = next(handle, None)
        if header is None:
            return rules
        for line in handle:
            fields = [field.strip() for field in line.rstrip("\n").split(",", 1)]
            if len(fields) != 2:
                continue
            lhs = [clean_token(token) for token in fields[0].split() if clean_token(token)]
            rhs = [clean_token(token) for token in fields[1].split() if clean_token(token)]
            if lhs:
                rules.append((lhs, rhs))
    return rules
