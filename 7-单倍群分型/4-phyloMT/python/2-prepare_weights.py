"""从 phylotree_index_withacc.json 生成 weights.txt

用法：
    python3 python/2-prepare_weights.py

输入：data/trees/phylotree-new-rcrs/17.2/phylotree_index_withacc.json
输出：tmp/weights.txt

权重公式（复现 Haplogrep 官方逻辑）：
    weight = round(SCALE * log(max_count / count) + 1, 1)
    其中 SCALE = 1.71，max_count 为出现频次最高的突变的频次。

Token 格式：
    SNP  T152C  → 152C  （去掉参考碱基）
    back T152C! → 152C! （去掉参考碱基，保留单个 !）
    del/ins 保持原样（小写 d），如 8281d、573.XC

rCRS和RSRS可以共用同一个 weights.txt。
所有单倍群的突变定义完全相同——突变标注全部使用rCRS 坐标系。
两者的区别不在于突变的表示方式，而在于：
  - rcrs 目录用 rCRS 作为比对参考（fasta 为 rCRS）
  - rsrs 目录额外提供了 rsrs.fasta 供坐标转换辅助

对于 RSRS→rCRS 差异位点（正向突变不在树节点中显式出现）：
    count = back_mutation_count + 1
"""

from __future__ import annotations

import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
JSON_PATH = PROJECT_DIR / "data" / "trees" / "phylotree-new-rcrs" / "17.3" / "phylotree_index_withacc.json"
OUT_PATH = PROJECT_DIR / "tmp" / "weights.txt"

SCALE = 1.71


def parse_token(raw: str) -> tuple[str, int]:
    """解析单个突变字符串，返回 (token, bang_count)。

    bang_count >= 1 表示 back mutation，统一折叠为单个 !。
    括号内的突变（不确定）也参与统计。
    """
    m = raw.strip("()")
    bangs = m.count("!")
    base = m.rstrip("!").upper()

    snp = re.match(r"^[ACGT](\d+[ACGT])$", base)
    if snp:
        token = snp.group(1)
    else:
        # del: 105-110D → 105-110d; XXXd → XXXd; ins/XC 保持
        token = re.sub(r"D$", "d", base)
        token = re.sub(r"-(\d+)D$", r"-\1d", token)

    if bangs >= 1:
        token = token + "!"
    return token, bangs


def count_mutations(haplogroups: dict) -> tuple[Counter, Counter]:
    """统计正向突变和 back mutation 各自的节点出现次数。"""
    forward: Counter = Counter()
    back: Counter = Counter()

    for node in haplogroups.values():
        muts = node.get("mutations", "").strip()
        if not muts:
            continue
        for raw in muts.split():
            token, bangs = parse_token(raw)
            if bangs >= 1:
                # token 已包含 !，去掉 ! 得到基础 token
                base_tok = token.rstrip("!")
                back[base_tok] += 1
            else:
                forward[token] += 1

    return forward, back


def compute_count(token: str, forward: Counter, back: Counter) -> int:
    """计算突变的有效频次。

    有效频次 = 正向出现次数 + back mutation 次数。
    对于 RSRS→rCRS 差异位点（正向出现次数为 0），额外 +1 代表树根处的隐式正向突变。
    """
    f = forward[token]
    b = back[token]
    if f == 0 and b == 0:
        return 1
    return f + b + (1 if f == 0 else 0)


def weight_from_count(count: int, max_count: int) -> float:
    return round(SCALE * math.log(max_count / count) + 1, 1)


def collect_all_tokens(forward: Counter, back: Counter) -> set[str]:
    """收集所有需要写入的基础 token（正向突变 token）。"""
    tokens: set[str] = set()
    for tok in forward:
        tokens.add(tok)
    for tok in back:
        # back mutation 的基础 token 对应正向突变
        tokens.add(tok)
        # back mutation 本身也需要写入（带 !）
        tokens.add(tok + "!")
    return tokens


def main() -> None:
    if not JSON_PATH.exists():
        print(f"错误：文件不存在: {JSON_PATH}", file=sys.stderr)
        sys.exit(1)

    print(f"读取 JSON: {JSON_PATH}")
    with JSON_PATH.open("r", encoding="utf-8") as f:
        data = json.load(f)

    haplogroups: dict = data["haplogroups"]
    print(f"节点数: {len(haplogroups)}")

    forward, back = count_mutations(haplogroups)
    print(f"正向突变 token 数: {len(forward)}")
    print(f"back mutation 基础 token 数: {len(back)}")

    # 收集所有 token（正向 + back mutation 自身）
    all_base_tokens = collect_all_tokens(forward, back)

    # 计算每个 token 的有效频次
    token_counts: dict[str, int] = {}
    for tok in all_base_tokens:
        if tok.endswith("!"):
            # back mutation token 本身
            base = tok.rstrip("!")
            cnt = back[base]
            if cnt > 0:
                token_counts[tok] = cnt
        else:
            cnt = compute_count(tok, forward, back)
            token_counts[tok] = cnt

    # max_count（正向突变中最高频次，用于权重归一化）
    max_count = max(token_counts[t] for t in token_counts if not t.endswith("!"))
    best_tok = max(
        ((t, c) for t, c in token_counts.items() if not t.endswith("!")),
        key=lambda x: x[1],
    )[0]
    print(f"最高频次: {max_count}（对应 token: {best_tok}）")

    # 按频次降序排列
    sorted_tokens = sorted(token_counts.items(), key=lambda x: (-x[1], x[0]))

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", encoding="utf-8") as f:
        for tok, cnt in sorted_tokens:
            w = weight_from_count(cnt, max_count)
            f.write(f"{tok}\t{w}\t{float(cnt)}\t{w}\t{w}\t{w}\n")

    print(f"写入完成: {OUT_PATH}（{len(sorted_tokens)} 行）")


if __name__ == "__main__":
    main()
