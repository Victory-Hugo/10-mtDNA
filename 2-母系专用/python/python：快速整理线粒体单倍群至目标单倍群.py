#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 如下脚本可以在普通电脑10秒内查询五万个单倍群

from __future__ import annotations
from pathlib import Path
from functools import lru_cache
import re
import sys
import argparse
import pandas as pd

# =========================== 通用函数 =========================== #
_pat_strip = re.compile(r'[\*\-]|(\+.*)')

def strip_suffix(hap: str) -> str:
    """移除 *、- 以及 + 之后的修饰位点"""
    return _pat_strip.sub('', hap)

def validate_input(df: pd.DataFrame, required: list[str]) -> None:
    if not set(required) <= set(df.columns):
        raise ValueError(f"输入缺少必要列: {required}")

# =========================== 核心步骤 =========================== #
def apply_primary_correction(id_hap_path: Path,
                             correction_path: Path,
                             save_path: Path) -> pd.DataFrame:
    id_hap = pd.read_csv(id_hap_path, sep='\t')
    validate_input(id_hap, ["ID", "Haplogroup"])

    corr_map = pd.read_csv(correction_path).set_index("Origin")["Revised"]
    id_hap["Haplogroup"] = id_hap["Haplogroup"].replace(corr_map.to_dict())

    save_path.parent.mkdir(parents=True, exist_ok=True)
    id_hap.to_csv(save_path, sep='\t', index=False)
    return id_hap

def build_parent_map(tree_path: Path) -> dict[str, str | None]:
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
                        forward_file: Path,
                        reverse_file: Path,
                        not_found_file: Path) -> None:
    """
    生成 Level_0..n / Level_n..0 两个表 + 未找到列表
    只输出包含数据的列，不输出全空列
    """
    id_hap = id_hap.copy()
    id_hap["lineage"] = id_hap["Haplogroup"].map(lineage_func)

    # 为每个样本创建正序和逆序的行
    forward_rows = []
    reverse_rows = []
    
    for _, row in id_hap.iterrows():
        lineage = row["lineage"]
        fwd = [row["ID"]] + list(lineage)
        rev = [row["ID"]] + list(reversed(lineage))
        
        forward_rows.append(fwd)
        reverse_rows.append(rev)
    
    # 找到最长的行
    max_fwd_len = max(len(row) for row in forward_rows)
    max_rev_len = max(len(row) for row in reverse_rows)
    
    # 生成列名
    level_cols_fwd = ["ID"] + [f"Level_{i}" for i in range(max_fwd_len - 1)]
    level_cols_rev = ["ID"] + [f"Level_{i}" for i in range(max_rev_len - 1)]
    
    # 创建 DataFrame 并只保留需要的列数
    df_fwd = pd.DataFrame(forward_rows, columns=level_cols_fwd[:max_fwd_len])
    df_rev = pd.DataFrame(reverse_rows, columns=level_cols_rev[:max_rev_len])

    forward_file.parent.mkdir(parents=True, exist_ok=True)
    reverse_file.parent.mkdir(parents=True, exist_ok=True)
    not_found_file.parent.mkdir(parents=True, exist_ok=True)

    df_fwd.to_csv(forward_file,  sep='\t', index=False)
    df_rev.to_csv(reverse_file,  sep='\t', index=False)

    # 未找到：lineage 为空
    missed = id_hap[id_hap["lineage"].str.len() == 0][["ID", "Haplogroup"]]
    if not missed.empty:
        missed.to_csv(not_found_file, sep='\t', index=False)

def match_to_target(reverse_file: Path,
                    target_path: Path,
                    rough_fix_path: Path,
                    out_final: Path,
                    out_unmatched: Path) -> None:
    df_rev = pd.read_csv(reverse_file, sep='\t')
    validate_input(df_rev, ["ID"])

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

    # 输出结果
    out_final.parent.mkdir(parents=True, exist_ok=True)
    out_unmatched.parent.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(matched, columns=["ID", "Haplogroup_PCA"]).to_csv(out_final, sep='\t', index=False)
    pd.DataFrame(unmatched, columns=["ID", "Haplogroup"]).to_csv(out_unmatched, sep='\t', index=False)

    # 特殊 L* 处理
    if unmatched and all(h.startswith("L") for _, h in unmatched if h):
        df_un = pd.read_csv(out_unmatched, sep='\t')
        df_un["Haplogroup_PCA"] = df_un["Haplogroup"].str[:2]
        pd.concat([pd.read_csv(out_final, sep='\t'), df_un[["ID", "Haplogroup_PCA"]]]
                  ).to_csv(out_final, sep='\t', index=False)

# =========================== 参数解析与入口 =========================== #
def parse_args():
    parser = argparse.ArgumentParser(
        description="mtDNA 单倍群管线：纠正名称 -> 构建谱系 -> 匹配目标",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 目录级别（若单个文件未指定，则按默认文件名基于目录推导）
    parser.add_argument("--output-dir", type=Path,
        default=Path("/mnt/f/3_陈静项目/2_单倍群分析/output"),
        help="输出目录（决定默认输出文件路径与默认的 ID_Hap.txt 位置）")
    parser.add_argument("--src-dir", type=Path,
        default=Path("/mnt/f/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用"),
        help="源码/资源目录（决定默认纠正表、phylotree、目标、粗纠正文件的路径）")

    # 输入与资源文件（可覆盖）
    parser.add_argument("--file-id-hap", type=Path, default=None, help="原始 ID-Hap（两列，制表符）")
    parser.add_argument("--file-correction", type=Path, default=None, help="Haplogrep 单倍群分型订正表.csv")
    parser.add_argument("--file-phylotree", type=Path, default=None, help="线粒体单倍群 phylotree（以制表符缩进表示层级）")
    parser.add_argument("--file-target", type=Path, default=None, help="目标列表（单列）")
    parser.add_argument("--file-rough-fix", type=Path, default=None, help="粗纠正表（两列，制表符）")

    # 输出文件（可覆盖）
    parser.add_argument("--out-fixed-hap", type=Path, default=None, help="订正之后的单倍群名称.txt")
    parser.add_argument("--out-forward-level", type=Path, default=None, help="正序等级.txt")
    parser.add_argument("--out-reversed-level", type=Path, default=None, help="逆序等级.txt")
    parser.add_argument("--out-not-found", type=Path, default=None, help="没有查询到请核实.txt")
    parser.add_argument("--out-final", type=Path, default=None, help="最终.txt")
    parser.add_argument("--out-unmatched", type=Path, default=None, help="无法处理.txt")

    return parser.parse_args()

def resolve_paths(args):
    """用目录默认值补齐未显式指定的文件路径"""
    OUTPUT_DIR = args.output_dir
    SRC_DIR    = args.src_dir

    FILE_ID_HAP     = args.file_id_hap     or (OUTPUT_DIR / "ID_Hap.txt")
    FILE_CORRECTION = args.file_correction or (SRC_DIR / "Haplogrep单倍群分型订正表.csv")
    FILE_PHYLOTREE  = args.file_phylotree  or (SRC_DIR / "线粒体单倍群phylotree(version17)2025年3月12日.txt")
    FILE_TARGET     = args.file_target     or (SRC_DIR / "目标_YuChunLi.txt")
    FILE_ROUGH_FIX  = args.file_rough_fix  or (SRC_DIR / "错误纠正_YuChunLi.txt")

    OUT_FIXED_HAP      = args.out_fixed_hap      or (OUTPUT_DIR / "订正之后的单倍群名称.txt")
    OUT_FORWARD_LEVEL  = args.out_forward_level  or (OUTPUT_DIR / "正序等级.txt")
    OUT_REVERSED_LEVEL = args.out_reversed_level or (OUTPUT_DIR / "逆序等级.txt")
    OUT_NOT_FOUND      = args.out_not_found      or (OUTPUT_DIR / "没有查询到请核实.txt")
    OUT_FINAL          = args.out_final          or (OUTPUT_DIR / "最终.txt")
    OUT_UNMATCHED      = args.out_unmatched      or (OUTPUT_DIR / "无法处理.txt")

    return (FILE_ID_HAP, FILE_CORRECTION, FILE_PHYLOTREE, FILE_TARGET, FILE_ROUGH_FIX,
            OUT_FIXED_HAP, OUT_FORWARD_LEVEL, OUT_REVERSED_LEVEL, OUT_NOT_FOUND, OUT_FINAL, OUT_UNMATCHED)

# =========================== 主调用 =========================== #
if __name__ == '__main__':
    try:
        args = parse_args()
        (FILE_ID_HAP, FILE_CORRECTION, FILE_PHYLOTREE, FILE_TARGET, FILE_ROUGH_FIX,
         OUT_FIXED_HAP, OUT_FORWARD_LEVEL, OUT_REVERSED_LEVEL, OUT_NOT_FOUND, OUT_FINAL, OUT_UNMATCHED) = resolve_paths(args)

        # 1) 基础纠正
        df_hap = apply_primary_correction(FILE_ID_HAP, FILE_CORRECTION, OUT_FIXED_HAP)

        # 2) lineage 表
        lineage = build_lineage_functions(build_parent_map(FILE_PHYLOTREE))
        output_level_tables(df_hap, lineage, OUT_FORWARD_LEVEL, OUT_REVERSED_LEVEL, OUT_NOT_FOUND)

        # 3) 目标匹配
        match_to_target(OUT_REVERSED_LEVEL, FILE_TARGET, FILE_ROUGH_FIX, OUT_FINAL, OUT_UNMATCHED)

        print("全部流程运行完毕 ✅")
    except Exception as e:
        print(
            "💥 运行中断: 一个可能的原因是生成的 '无法处理.txt' 文档中存在除了 L 单倍群之外的单倍群；"
            "请手动在 Excel 等软件内处理这些特殊样本。\n详细错误: ", e,
            file=sys.stderr
        )
        sys.exit(1)
