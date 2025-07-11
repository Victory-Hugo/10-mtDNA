#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#! 如下脚本可以在普通电脑10秒内查询五万个单倍群

from __future__ import annotations
from pathlib import Path
from functools import lru_cache
import re
import sys
import pandas as pd

# =========================== 1. 统一配置 =========================== #
OUTPUT_DIR   = Path("/mnt/f/3_陈静项目/2_单倍群分析/output")
SRC_DIR     = Path("/mnt/f/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用")

FILE_ID_HAP          = OUTPUT_DIR / "ID_Hap.txt"      # 原始 ID‑Hap（两列）
FILE_CORRECTION      = SRC_DIR / "Haplogrep单倍群分型订正表.csv"
FILE_PHYLOTREE       = SRC_DIR / "线粒体单倍群phylotree(version17)2025年3月12日.txt"
FILE_TARGET          = SRC_DIR / "目标_YuChunLi.txt"
FILE_ROUGH_FIX       = SRC_DIR / "错误纠正_YuChunLi.txt"

OUT_FIXED_HAP        = OUTPUT_DIR / "订正之后的单倍群名称.txt"
OUT_FORWARD_LEVEL    = OUTPUT_DIR / "正序等级.txt"
OUT_REVERSED_LEVEL   = OUTPUT_DIR / "逆序等级.txt"
OUT_NOT_FOUND        = OUTPUT_DIR / "没有查询到请核实.txt"
OUT_FINAL            = OUTPUT_DIR / "最终.txt"
OUT_UNMATCHED        = OUTPUT_DIR / "无法处理.txt"

# =========================== 2. 通用函数 =========================== #
_pat_strip = re.compile(r'[\*\-]|(\+.*)')

def strip_suffix(hap: str) -> str:
    """移除 *、- 以及 + 之后的修饰位点"""
    return _pat_strip.sub('', hap)

def validate_input(df: pd.DataFrame, required: list[str]) -> None:
    if not set(required) <= set(df.columns):
        raise ValueError(f"输入缺少必要列: {required}")

# =========================== 3. 第一步：纠正 haplogroup 名称 =========================== #
def apply_primary_correction(id_hap_path=FILE_ID_HAP,
                             correction_path=FILE_CORRECTION,
                             save_path=OUT_FIXED_HAP) -> pd.DataFrame:
    id_hap = pd.read_csv(id_hap_path, sep='\t')
    validate_input(id_hap, ["ID", "Haplogroup"])

    corr_map = pd.read_csv(correction_path).set_index("Origin")["Revised"]
    id_hap["Haplogroup"] = id_hap["Haplogroup"].replace(corr_map.to_dict())

    id_hap.to_csv(save_path, sep='\t', index=False)
    return id_hap

# =========================== 4. 第二步：构建 PhyloTree 追溯链 =========================== #
def build_parent_map(tree_path=FILE_PHYLOTREE) -> dict[str, str | None]:
    """
    将 phylotree 结构(以制表符缩进表示层级)解析成 child->parent 字典
    """
    parent_map: dict[str, str | None] = {}
    stack: list[tuple[int, str]] = []                   # [(level, hap)]
    with open(tree_path, encoding='utf-8') as f:
        for raw in f:
            if not raw.strip():
                continue
            level = raw.count('\t')
            hap   = raw.strip()
            while stack and stack[-1][0] >= level:
                stack.pop()
            parent = stack[-1][1] if stack else None
            parent_map[hap] = parent
            stack.append((level, hap))
    return parent_map

def build_lineage_functions(parent_map):
    @lru_cache(maxsize=None)
    def lineage(hap: str) -> tuple[str, ...]:
        """返回从根到该 hap 的路径(含自身)"""
        path = []
        while hap:
            path.append(hap)
            hap = parent_map.get(hap)
        return tuple(reversed(path))

    return lineage

def output_level_tables(id_hap: pd.DataFrame,
                        lineage_func,
                        forward_file=OUT_FORWARD_LEVEL,
                        reverse_file=OUT_REVERSED_LEVEL,
                        not_found_file=OUT_NOT_FOUND) -> None:
    """
    生成 Level_0..n / Level_n..0 两个表 + 未找到列表
    """
    # —— 计算 lineage 并统计最长深度 ——
    id_hap["lineage"] = id_hap["Haplogroup"].map(lineage_func)
    max_len = id_hap["lineage"].str.len().max()

    # —— 制作正序 & 逆序列 —— 
    level_cols_fwd = [f"Level_{i}"        for i in range(max_len)]
    level_cols_rev = [f"Level_{i}"        for i in range(max_len-1, -1, -1)]

    def pad(seq: tuple[str, ...]) -> list[str]:
        seq = list(seq)
        return seq + [""] * (max_len - len(seq))

    id_hap[level_cols_fwd] = id_hap["lineage"].apply(lambda x: pd.Series(pad(x)))
    id_hap[level_cols_rev] = id_hap["lineage"].apply(lambda x: pd.Series(pad(tuple(reversed(x)))))

    # —— 保存 —— 
    id_hap[["ID"] + level_cols_fwd].to_csv(forward_file,  sep='\t', index=False)
    id_hap[["ID"] + level_cols_rev].to_csv(reverse_file,  sep='\t', index=False)

    # —— 未找到：lineage 为空 —— 
    missed = id_hap[id_hap["lineage"].str.len() == 0][["ID", "Haplogroup"]]
    if not missed.empty:
        missed.to_csv(not_found_file, sep='\t', index=False)

# =========================== 5. 第三步：匹配目标列表 =========================== #
def match_to_target(reverse_file=OUT_REVERSED_LEVEL,
                    target_path=FILE_TARGET,
                    rough_fix_path=FILE_ROUGH_FIX,
                    out_final=OUT_FINAL,
                    out_unmatched=OUT_UNMATCHED) -> None:
    df_rev = pd.read_csv(reverse_file, sep='\t')
    validate_input(df_rev, ["ID"])

    # 目标集合 & 粗纠正 map
    target_set = set(pd.read_csv(target_path, header=None, names=["Hap"])["Hap"].map(strip_suffix))
    rough_map  = pd.read_csv(rough_fix_path, sep='\t', header=None,
                             names=['Original', 'Correction']).set_index("Original")["Correction"]

    level_cols = df_rev.columns[1:]          # ID 后所有列

    matched, unmatched = [], []
    for _, row in df_rev.iterrows():
        seq = [strip_suffix(str(h)) for h in row[level_cols] if pd.notna(h) and str(h).strip()]
        hit = next((h for h in seq if h in target_set), None)
        if hit:
            matched.append((row["ID"], rough_map.get(hit, hit)))
        else:
            fall_back = seq[0] if seq else ""
            unmatched.append((row["ID"], fall_back))

    # —— 输出结果 —— 
    pd.DataFrame(matched, columns=["ID", "Haplogroup_PCA"]).to_csv(out_final, sep='\t', index=False)
    pd.DataFrame(unmatched, columns=["ID", "Haplogroup"]).to_csv(out_unmatched, sep='\t', index=False)

    # —— 特殊 L* 处理 —— 
    if unmatched and all(h.startswith("L") for _, h in unmatched if h):
        df_un = pd.read_csv(out_unmatched, sep='\t')
        df_un["Haplogroup_PCA"] = df_un["Haplogroup"].str[:2]
        pd.concat([pd.read_csv(out_final, sep='\t'), df_un[["ID", "Haplogroup_PCA"]]]
                  ).to_csv(out_final, sep='\t', index=False)

# =========================== 6. 主调用 =========================== #
if __name__ == '__main__':
    try:
        # 1) 基础纠正
        df_hap = apply_primary_correction()

        # 2) lineage 表
        lineage = build_lineage_functions(build_parent_map())
        output_level_tables(df_hap, lineage)

        # 3) 目标匹配
        match_to_target()

        print("全部流程运行完毕 ✅")
    except Exception as e:
        print(f"💥 运行中断: 一个可能的原因是生成的 '无法处理.txt' 文档中存在除了 L 单倍群之外的单倍群\n请手动在Excel等软件内处理这些特殊的样本", e, file=sys.stderr)
        sys.exit(1)
