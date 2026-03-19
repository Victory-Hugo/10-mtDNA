#!/usr/bin/env python3
"""
把 txt 文件中含下划线的新单倍群合并到 phylotree_build17.csv。

新单倍群格式："基础单倍群_突变1_突变2..."
- haplogroup = 原始全名（含下划线）
- parent     = 在 txt 树中的直接父节点
- mutations  = 下划线分割后的部分（第2个起，空格连接）
- level      = 父节点的 level + 1

输出 CSV 按 txt 文件的树顺序排列（既保留原CSV顺序又插入新节点）。
"""

import csv
import re

TXT  = "线粒体单倍群phylotree(version17)2025年3月12日.txt"
CSV  = "phylotree_build17.csv"
OUT  = "phylotree_build17.csv"


# ── 1. 读取 txt，建立树结构 ──────────────────────────────────────────────
print("读取 txt 树结构...")
txt_nodes = []   # [(indent_level, name)]
with open(TXT, encoding='utf-8') as f:
    for raw in f:
        line = raw.rstrip('\n')
        name = line.lstrip('\t')
        indent = len(line) - len(name)
        name = name.strip()
        if name:
            txt_nodes.append((indent, name))

# 用栈推导 parent
parent_of = {}   # name -> parent_name（用 txt 顺序中的直接上级）
stack = []       # [(indent, name)]
for indent, name in txt_nodes:
    while stack and stack[-1][0] >= indent:
        stack.pop()
    parent_of[name] = stack[-1][1] if stack else ''
    stack.append((indent, name))

print(f"  txt 节点总数：{len(txt_nodes)}")


# ── 2. 读取已有 CSV ───────────────────────────────────────────────────────
print("读取已有 CSV...")
existing_rows = []
csv_level = {}   # name -> level（int）
with open(CSV, encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        existing_rows.append(row)
        csv_level[row['haplogroup']] = int(row['level'])

csv_names = set(csv_level.keys())
print(f"  CSV 节点总数：{len(existing_rows)}")


# ── 3. 计算新节点的 level ─────────────────────────────────────────────────
# 若父节点在 CSV 中 → level = csv_level[parent] + 1
# 若父节点也是新节点 → 递归计算
# txt 节点的"缩进层" 用作兜底（当 parent 信息不可追溯时）
txt_indent = {name: indent for indent, name in txt_nodes}

_memo = {}
def compute_level(name):
    if name in _memo:
        return _memo[name]
    if name in csv_level:
        _memo[name] = csv_level[name]
        return csv_level[name]
    par = parent_of.get(name, '')
    if not par:
        result = 0
    else:
        result = compute_level(par) + 1
    _memo[name] = result
    return result


# ── 4. 找出所有含下划线的新单倍群 ─────────────────────────────────────────
new_haplos = [name for _, name in txt_nodes if '_' in name and name not in csv_names]
print(f"  含下划线新单倍群：{len(new_haplos)}")

# 构建新节点的行字典
new_row_of = {}
for name in new_haplos:
    parts = name.split('_')
    mutations = ' '.join(parts[1:])  # 下划线后的所有部分
    parent = parent_of.get(name, '')
    level  = compute_level(name)
    new_row_of[name] = {
        'haplogroup': name,
        'level':      level,
        'parent':     parent,
        'mutations':  mutations,
    }

# 调试：检查几个样例
samples = list(new_haplos)[:5]
print("\n新节点样例：")
for n in samples:
    r = new_row_of[n]
    print(f"  {r['haplogroup']:40s}  level={r['level']:2d}  parent={r['parent']:30s}  muts={r['mutations']}")


# ── 5. 按 txt 树顺序合并输出 ──────────────────────────────────────────────
# 遍历 txt_nodes 顺序：
#   - 若是已有 CSV 节点 → 输出原来的行
#   - 若是新下划线节点 → 输出新行
#   - 其它（如 MT-MRCA、M7a2a4 等不在 CSV 且无下划线）→ 跳过

print("\n合并输出...")

# 用字典快速查找已有 CSV 行
csv_row_of = {r['haplogroup']: r for r in existing_rows}

output_rows = []
seen = set()

for _, name in txt_nodes:
    if name in seen:
        continue
    seen.add(name)

    if name in csv_row_of:
        output_rows.append(csv_row_of[name])
    elif name in new_row_of:
        output_rows.append(new_row_of[name])
    # 其它不在 CSV、也不是下划线新节点 → 跳过

# 检查是否有 CSV 中的节点没被 txt 覆盖（理论上应该很少）
covered = {r['haplogroup'] for r in output_rows}
not_covered = [name for name in csv_names if name not in covered]
if not_covered:
    print(f"\n  警告：{len(not_covered)} 个 CSV 节点未出现在 txt 树中，追加到末尾：")
    for n in not_covered[:10]:
        print(f"    {n}")
    for name in not_covered:
        output_rows.append(csv_row_of[name])

print(f"  输出行数：{len(output_rows)}（原 {len(existing_rows)} + 新 {len(new_haplos)}）")


# ── 6. 写入 CSV ───────────────────────────────────────────────────────────
FIELDS = ['haplogroup', 'level', 'parent', 'mutations']
with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(output_rows)

print(f"\n完成！已写入 {OUT}")
