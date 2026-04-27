"""mtDNA 变异表示与规则处理。"""

from __future__ import annotations

import re
from collections import Counter


BACK_MUTATION_SUFFIX = "!"


_DEL_RE = re.compile(r"^(\d+(?:\.\d+)?)[Dd][Ee][Ll]$")


def clean_token(token: str) -> str:
    """清理单个变异字符串，统一删除后缀为 d。"""

    token = token.strip().strip('"').strip()
    m = _DEL_RE.match(token)
    if m:
        token = m.group(1) + "d"
    return token


def is_back_mutation(token: str) -> bool:
    """判断是否为返祖事件。"""

    return clean_token(token).endswith(BACK_MUTATION_SUFFIX)


def base_token(token: str) -> str:
    """去掉返祖标记后的基本形式。"""

    token = clean_token(token)
    if token.endswith(BACK_MUTATION_SUFFIX):
        return token[:-1]
    return token


def mutation_key(token: str) -> str:
    """提取用于覆盖同位点状态的键。"""

    token = base_token(token)
    if "." in token:
        anchor, rest = token.split(".", 1)
        match = re.match(r"([0-9X]+)(.*)", rest)
        index = match.group(1) if match else rest
        return f"ins:{anchor}:{index}"
    match = re.match(r"(\d+)(.*)", token)
    if not match:
        return f"raw:{token}"
    position = match.group(1)
    suffix = match.group(2)
    if suffix == "d":
        return f"del:{position}"
    return f"sub:{position}"


def mutation_position(token: str) -> int | None:
    """提取变异位点。"""

    token = base_token(token)
    match = re.match(r"(\d+)", token)
    if not match:
        return None
    return int(match.group(1))


def apply_variants_to_profile(
    parent_variants: list[str],
    own_variants: list[str],
) -> list[str]:
    """沿路径累计树节点变异。"""

    profile: dict[str, str] = {mutation_key(token): token for token in parent_variants}
    for raw_token in own_variants:
        token = clean_token(raw_token)
        if not token:
            continue
        key = mutation_key(token)
        if is_back_mutation(token):
            profile.pop(key, None)
        else:
            profile[key] = base_token(token)
    return list(profile.values())


def weight_of(token: str, weights: dict[str, float]) -> float:
    """获取单个变异权重。"""

    token = clean_token(token)
    if token in weights:
        return weights[token]
    token = base_token(token)
    return weights.get(token, 0.0)


def sum_weights(tokens: list[str], weights: dict[str, float]) -> float:
    """累计变异权重。"""

    return round(sum(weight_of(token, weights) for token in tokens), 10)


def ordered_unique(tokens: list[str]) -> list[str]:
    """保持输入顺序去重。"""

    seen: set[str] = set()
    result: list[str] = []
    for token in tokens:
        token = clean_token(token)
        if token and token not in seen:
            seen.add(token)
            result.append(token)
    return result


_MT_LENGTH = 16569


def parse_range_to_positions(range_text: str) -> "set[int] | None":
    """将范围文本解析为覆盖位点集合，全基因组范围返回 None。

    支持：
    - X-Y  普通范围（X <= Y）
    - X-Y  环形范围（X > Y，即 X-16569 加 1-Y）
    - X    单个位点
    分隔符支持空格和分号。
    """

    positions: set[int] = set()
    for seg in re.split(r"[\s;]+", range_text.strip()):
        seg = seg.strip()
        if not seg:
            continue
        m = re.match(r"^(\d+)-(\d+)$", seg)
        if m:
            start, end = int(m.group(1)), int(m.group(2))
            if start <= end:
                positions.update(range(start, end + 1))
            else:
                positions.update(range(start, _MT_LENGTH + 1))
                positions.update(range(1, end + 1))
        elif re.match(r"^\d+$", seg):
            positions.add(int(seg))
    if not positions or positions >= set(range(1, _MT_LENGTH + 1)):
        return None
    return positions


def position_covered(token: str, covered: "set[int] | None") -> bool:
    """判断变异位点是否在覆盖范围内。"""

    if covered is None:
        return True
    pos = mutation_position(token)
    return pos is not None and pos in covered


def apply_alignment_rules(
    variants: list[str],
    rules: list[tuple[list[str], list[str]]],
) -> list[str]:
    """按规则表重写 FASTA 对齐结果。"""

    current = ordered_unique(variants)
    changed = True
    while changed:
        changed = False
        for lhs, rhs in rules:
            counter = Counter(current)
            if all(counter[token] >= lhs.count(token) for token in lhs):
                next_tokens = list(current)
                for token in lhs:
                    next_tokens.remove(token)
                next_tokens.extend(rhs)
                current = ordered_unique(next_tokens)
                changed = True
                break
    return current
