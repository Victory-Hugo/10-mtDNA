"""为 phylotree-new 树目录准备运行所需的全部辅助文件。

用法：
    python3 python/prepare_new_tree.py

策略：
  新版 JSON 树使用 RSRS 坐标系定义变异，但输入数据（HSD/VCF）
  使用 rCRS 坐标系。因此：
  - fasta 使用 rCRS（用于 FASTA 比对及 HSD/VCF 变异提取）
  - rsrsFasta 指向 RSRS fasta（供 load_tree_bundle 在加载时做坐标转换）
  - weights 使用 rCRS@17.2 的权重文件（与 rCRS 坐标系匹配）
  - rules/gff/aac 从 rCRS@17.2 复制

执行后生成 data/trees/phylotree-new/1.0/ 目录。
"""

from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# 路径常量
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]
RCRS_DIR = PROJECT_DIR / "data" / "trees" / "phylotree-rcrs" / "17.2"
RSRS_DIR = PROJECT_DIR / "data" / "trees" / "phylotree-rsrs" / "17.1"
JSON_SRC = PROJECT_DIR / "data" / "trees" / "phylotree-new" / "phylotree_index_withacc.json"
DST_DIR = PROJECT_DIR / "data" / "trees" / "phylotree-new" / "1.0"

# 从 rCRS 17.2 复制的文件（用于比对和评分）
RCRS_FILES = [
    "rcrs.fasta",
    "rcrs.fasta.fai",
    "rcrs.fasta.amb",
    "rcrs.fasta.ann",
    "rcrs.fasta.bwt",
    "rcrs.fasta.pac",
    "rcrs.fasta.sa",
    "rcrs.dict",
    "rules.csv",
    "mtdna.gff",
    "mtdna.aac.txt",
    "weights.txt",
]

# 从 RSRS 17.1 复制（仅供 RSRS→rCRS 坐标转换使用）
RSRS_FILES = [
    "rsrs.fasta",
]

TREE_YAML = """\
id: phylotree-new
version: 1.0
name: PhyloTree New (rCRS)
category: rCRS (Human mtDNA)
description: New phylotree JSON (RSRS-based internally, converted to rCRS coordinates at load time)
tree: phylotree_index_withacc.json
weights: weights.txt
fasta: rcrs.fasta
rsrsFasta: rsrs.fasta
alignmentRules: rules.csv
"""

# ---------------------------------------------------------------------------
# 变异格式规范化（与 tree.py 中保持一致）
# ---------------------------------------------------------------------------
_RANGE_DEL_RE = re.compile(r"^(\d+)-(\d+)d$", re.IGNORECASE)


def normalize_json_token(raw_token: str) -> list[str]:
    """将新版 JSON 变异 token 转为旧版权重文件兼容格式（RSRS 坐标系）。"""
    token = raw_token.strip()
    if not token or token == "reserved":
        return []
    if token.startswith("(") and token.endswith(")"):
        return []

    m = _RANGE_DEL_RE.match(token)
    if m:
        start, end = int(m.group(1)), int(m.group(2))
        return [f"{pos}d" for pos in range(start, end + 1)]

    if token[0].isalpha():
        rest = token[1:]
        if rest.endswith("!!"):
            body = rest[:-2]
            if body and body[-1].lower() in "acgt":
                body = body[:-1] + body[-1].upper()
            return [body + "!!"]
        if rest.endswith("!"):
            body = rest[:-1]
            if body and body[-1].lower() in "acgt":
                body = body[:-1] + body[-1].upper()
            return [body + "!"]
        if rest and rest[-1].lower() not in "acgt":
            return [rest]
        if rest and rest[-1].isalpha():
            return [rest[:-1] + rest[-1].upper()]
        return [rest]

    token = token.strip('"').strip()
    del_m = re.match(r"^(\d+(?:\.\d+)?)[Dd][Ee][Ll]$", token)
    if del_m:
        token = del_m.group(1) + "d"
    return [token] if token else []


def collect_new_tree_rsrs_variants(json_path: Path) -> set[str]:
    """提取新版 JSON 树所有 RSRS 坐标系规范化变异 token 集合。"""
    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    variants: set[str] = set()
    for info in data["haplogroups"].values():
        for raw in info.get("mutations", "").split():
            for tok in normalize_json_token(raw):
                if tok:
                    variants.add(tok)
    return variants


def load_weights(weights_path: Path) -> dict[str, str]:
    """读取权重文件，键为 token，值为完整行（含所有列）。"""
    weights: dict[str, str] = {}
    with weights_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if not fields:
                continue
            token = fields[0].strip().strip('"').strip()
            del_m = re.match(r"^(\d+(?:\.\d+)?)[Dd][Ee][Ll]$", token)
            if del_m:
                token = del_m.group(1) + "d"
            weights[token] = line
    return weights


def build_rcrs_variant_set_from_rsrs(json_path: Path, rsrs_seq: str, rcrs_seq: str) -> set[str]:
    """从新版 JSON 树的 RSRS 坐标变异推导 rCRS 坐标变异集合。

    仅推导 RSRS-rCRS 差异位点对应的 rCRS 变异（供权重文件中缺失条目统计）。
    """
    from collections import defaultdict

    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    haplogroups = data["haplogroups"]

    # 差异位点
    diff_map: dict[int, tuple[str, str]] = {}
    for i, (rb, sb) in enumerate(zip(rcrs_seq, rsrs_seq)):
        if rb != sb and rb != "N" and sb != "N":
            diff_map[i + 1] = (rb, sb)

    # 先构建完整树的 path_variants（简化版，不处理 own_variants 只看节点级）
    # 这里做一个简化：对每个节点的 own mutations，如果在差异位点，
    # 转换后加入集合
    rcrs_variants: set[str] = set()
    for pos, (rcrs_base, rsrs_base) in diff_map.items():
        # 位于差异位点的 rCRS 变异：<pos><非rCRS碱基>
        # 这些碱基可能是 rsrs_base 或其他
        # 对于 RSRS ancestral：rCRS变异 = <pos><rsrs_base>
        rcrs_variants.add(f"{pos}{rsrs_base}")
    return rcrs_variants


def generate_weights(
    json_path: Path,
    rcrs_weights_path: Path,
    rsrs_seq: str,
    rcrs_seq: str,
    dst_path: Path,
) -> None:
    """生成新版树的权重文件：以 rCRS@17.2 权重为基础，补充缺失条目。"""
    rcrs_weights = load_weights(rcrs_weights_path)

    # 计算转换后可能出现的 rCRS 变异
    rcrs_extra = build_rcrs_variant_set_from_rsrs(json_path, rsrs_seq, rcrs_seq)
    missing_from_rcrs = rcrs_extra - set(rcrs_weights.keys())

    print(f"  rCRS@17.2 权重条目数: {len(rcrs_weights)}")
    print(f"  转换后新增 rCRS 变异（可能缺失权重）: {len(missing_from_rcrs)}")

    with dst_path.open("w", encoding="utf-8") as f:
        for token, line in rcrs_weights.items():
            f.write(line + "\n")
        for token in sorted(missing_from_rcrs):
            f.write(f"{token}\t1.0\t1.0\t1.0\t1.0\t1.0\n")


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def main() -> None:
    if not JSON_SRC.exists():
        print(f"错误：新版 JSON 文件不存在: {JSON_SRC}", file=sys.stderr)
        sys.exit(1)
    if not RCRS_DIR.exists():
        print(f"错误：rCRS 17.2 源目录不存在: {RCRS_DIR}", file=sys.stderr)
        sys.exit(1)
    if not RSRS_DIR.exists():
        print(f"错误：RSRS 17.1 源目录不存在: {RSRS_DIR}", file=sys.stderr)
        sys.exit(1)

    DST_DIR.mkdir(parents=True, exist_ok=True)
    print(f"目标目录: {DST_DIR}")

    # 复制 rCRS 辅助文件
    print("从 rCRS 17.2 复制文件…")
    for fname in RCRS_FILES:
        src = RCRS_DIR / fname
        dst = DST_DIR / fname
        if not src.exists():
            print(f"  跳过（不存在）: {fname}")
            continue
        shutil.copy2(src, dst)
        print(f"  已复制: {fname}")

    # 复制 RSRS fasta（供坐标转换）
    print("从 RSRS 17.1 复制 rsrs.fasta…")
    for fname in RSRS_FILES:
        src = RSRS_DIR / fname
        dst = DST_DIR / fname
        if not src.exists():
            print(f"  跳过（不存在）: {fname}")
            continue
        shutil.copy2(src, dst)
        print(f"  已复制: {fname}")

    # 复制 JSON 树文件
    json_dst = DST_DIR / "phylotree_index_withacc.json"
    shutil.copy2(JSON_SRC, json_dst)
    print("  已复制: phylotree_index_withacc.json")

    # 读取参考序列
    from Bio import SeqIO
    rsrs_seq = str(next(SeqIO.parse(DST_DIR / "rsrs.fasta", "fasta")).seq).upper()
    rcrs_seq = str(next(SeqIO.parse(DST_DIR / "rcrs.fasta", "fasta")).seq).upper()

    # 生成 weights.txt（基于 rCRS@17.2）
    print("生成 weights.txt …")
    generate_weights(
        JSON_SRC,
        RCRS_DIR / "weights.txt",
        rsrs_seq,
        rcrs_seq,
        DST_DIR / "weights.txt",
    )

    # 写 tree.yaml
    (DST_DIR / "tree.yaml").write_text(TREE_YAML, encoding="utf-8")
    print("  已写入: tree.yaml")

    print("\n完成！新版树目录已就绪：", DST_DIR)


if __name__ == "__main__":
    main()
