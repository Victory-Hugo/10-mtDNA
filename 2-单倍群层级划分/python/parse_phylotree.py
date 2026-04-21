#!/usr/bin/env python3
"""
解析 PhyloTree Build 17 HTM 文件，输出含层级关系和定义突变的 TSV。

列含义：
  Haplogroup  - 单倍群名称
  Level       - 在树中的列号（代表层级深度）
  Parent      - 父节点单倍群名
  Mutations   - 定义该单倍群的突变列表（空格分隔）

使用：
  python3 parse_phylotree.py "mtDNA tree Build 17.htm" > output.tsv
"""

import re
import csv
import sys
from bs4 import BeautifulSoup

# 单倍群名格所用的 CSS 类（rowspan=2、非登录号格）
HAPLO_CLASSES = {
    'xl11317826',  # 粗边框，较高层级
    'xl11417826',  # 黄色背景，某些特殊节点
    'xl11517826',  # 最常见的普通节点
    'xl11617826',
    'xl11717826',
    'xl11817826',  # 黄色背景
    'xl11917826',
    'xl13417826',  # 顶层（L0 等）
}

# 登录号列的类，需排除
ACCESSION_CLASS = 'xl13117826'

def clean_mutation_text(text: str) -> str:
    """清理突变文本：去除 BOM、不换行空格、多余空白。"""
    text = text.replace('\ufeff', '').replace('\xa0', '').replace('\u00a0', '')
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def extract_mutations_from_row(cells_after_haplo) -> str:
    """从单倍群格之后的格中提取突变列表。"""
    parts = []
    for td in cells_after_haplo:
        cls_list = td.get('class', [])
        cls_str = ' '.join(cls_list) if isinstance(cls_list, list) else cls_list

        # 遇到登录号列就停止
        if ACCESSION_CLASS in cls_str:
            break

        # Build 17 中 mutation 单元格样式并不稳定；只要在 accession 列之前且有内容，
        # 就视为当前单倍群的定义突变。
        text = td.get_text(' ')
        text = clean_mutation_text(text)
        if text:
            parts.append(text)

    return ' '.join(parts)


def is_haplo_cell(td) -> bool:
    """判断一个 td 是否为单倍群名格。"""
    if td.get('rowspan') != '2':
        return False
    cls_list = td.get('class', [])
    cls_set = set(cls_list) if isinstance(cls_list, list) else {cls_list}
    return bool(cls_set & HAPLO_CLASSES)


def parse_phylotree(filepath: str):
    """解析 HTM 文件，返回单倍群信息列表。"""
    print(f"[1/3] 读取文件 {filepath} ...", file=sys.stderr)
    with open(filepath, 'r', encoding='windows-1252', errors='replace') as f:
        content = f.read()

    print("[2/3] 解析 HTML（大文件，请稍候）...", file=sys.stderr)
    soup = BeautifulSoup(content, 'lxml')
    table = soup.find('table')
    if not table:
        print("错误：找不到 <table>", file=sys.stderr)
        return []

    rows = table.find_all('tr')
    print(f"      共 {len(rows)} 行", file=sys.stderr)

    haplogroups = []
    # depth_stack[col_idx] = haplogroup_name  记录每个深度上最近的单倍群
    depth_stack = {}

    print("[3/3] 提取单倍群数据...", file=sys.stderr)

    for tr in rows:
        cells = list(tr.find_all('td'))
        if not cells:
            continue

        col = 0           # 当前列号（累计colspan）
        haplo_name = None
        haplo_col = None
        cells_after = []  # 单倍群格之后的格列表

        found = False
        for i, td in enumerate(cells):
            colspan = int(td.get('colspan', 1))

            if not found:
                if is_haplo_cell(td):
                    text = clean_mutation_text(td.get_text(' '))
                    if text:
                        haplo_name = text
                        haplo_col = col
                        cells_after = cells[i + 1:]
                        found = True
            col += colspan

        if not haplo_name:
            continue

        # 跳过根节点（mt-MRCA）
        if 'mt-MRCA' in haplo_name or 'MRCA' in haplo_name:
            depth_stack[haplo_col] = haplo_name
            continue

        depth = haplo_col

        # 找最近的上级节点
        parent = None
        for d in range(depth - 1, -1, -1):
            if d in depth_stack:
                parent = depth_stack[d]
                break

        # 更新 depth_stack：更新当前深度，清除更深的记录
        depth_stack[depth] = haplo_name
        for d in list(depth_stack.keys()):
            if d > depth:
                del depth_stack[d]

        mutations = extract_mutations_from_row(cells_after)

        haplogroups.append({
            'Haplogroup': haplo_name,
            'Level': depth,
            'Parent': parent if parent else '',
            'Mutations': mutations,
        })

    return haplogroups


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法：python3 parse_phylotree.py <HTM文件路径> [输出TSV路径]", file=sys.stderr)
        sys.exit(1)

    fp = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else None

    data = parse_phylotree(fp)
    print(f"\n共提取 {len(data)} 个单倍群节点", file=sys.stderr)

    fieldnames = ['Haplogroup', 'Level', 'Parent', 'Mutations']

    if out_path:
        with open(out_path, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(data)
        print(f"已写入：{out_path}", file=sys.stderr)
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(data)
