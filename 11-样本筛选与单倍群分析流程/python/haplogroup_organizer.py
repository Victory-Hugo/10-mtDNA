#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
单倍群整理完整流程：按照Jupyter笔记本中的步骤执行

规范单倍群的代码直接来自笔记本，只修改路径部分。
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List
from functools import lru_cache

import pandas as pd


# ============================================================================
# 步骤1-5：提取原始单倍群，生成ID_Hap.txt（来自笔记本单元格1-6）
# ============================================================================
def extract_raw_haplogroups(input_excel: str) -> pd.DataFrame:
    """
    从Excel中提取原始单倍群数据
    对应笔记本的单元格3-6
    """
    print("  读取质量控制sheet...")
    df_全测 = pd.read_excel(input_excel, sheet_name='质量控制')
    df_全测 = df_全测.loc[:, ['SampleID', 'Haplogroup']]
    
    print("  读取芯片单倍群sheet...")
    df_芯片 = pd.read_excel(input_excel, sheet_name='芯片单倍群')
    df_芯片 = df_芯片.loc[:, ['SampleID', 'Haplogroup']]
    
    # 合并
    df_ID_hap = pd.concat([df_全测, df_芯片])
    df_ID_hap = df_ID_hap.rename(columns={'SampleID': 'ID'})
    
    # 去除haplogroup列首尾空格,使保持格式一致
    df_ID_hap['Haplogroup'] = df_ID_hap['Haplogroup'].str.strip()
    
    return df_ID_hap


# ============================================================================
# 步骤6-8：YuChunLi流程（来自笔记本单元格8 - 规范单倍群）
# ============================================================================
def run_yuchunli_workflow(id_hap_path: str, conf_dir: str, temp_dir: str):
    """
    YuChunLi方案流程
    代码直接来自笔记本单元格8，只改路径部分
    """
    print("\n" + "="*60)
    print("🔵 YuChunLi 流程")
    print("="*60)
    
    # 配置路径（只改这部分）
    OUTPUT_DIR = Path(temp_dir)
    SRC_DIR = Path(conf_dir)
    
    FILE_ID_HAP = Path(id_hap_path)
    FILE_CORRECTION = SRC_DIR / "Haplogrep单倍群分型订正表.csv"
    FILE_PHYLOTREE = SRC_DIR / "线粒体单倍群phylotree(version17)2025年3月12日.txt"
    FILE_TARGET = SRC_DIR / "目标_YuChunLi.txt"
    FILE_ROUGH_FIX = SRC_DIR / "错误纠正_YuChunLi.txt"
    
    OUT_FIXED_HAP = OUTPUT_DIR / "订正之后的单倍群名称_YuChunLi.txt"
    OUT_FORWARD_LEVEL = OUTPUT_DIR / "正序等级_YuChunLi.txt"
    OUT_REVERSED_LEVEL = OUTPUT_DIR / "逆序等级_YuChunLi.txt"
    OUT_NOT_FOUND = OUTPUT_DIR / "没有查询到请核实_YuChunLi.txt"
    OUT_FINAL = OUTPUT_DIR / "最终_YuChunLi.txt"
    OUT_UNMATCHED = OUTPUT_DIR / "无法处理_YuChunLi.txt"
    
    # ==================== 以下代码直接来自笔记本，无任何修改 ====================
    def validate_input(df: pd.DataFrame, required: list[str]) -> None:
        if not set(required) <= set(df.columns):
            raise ValueError(f"输入缺少必要列: {required}")

    # 第一步：纠正 haplogroup 名称
    def apply_primary_correction(id_hap_path=FILE_ID_HAP,
                                 correction_path=FILE_CORRECTION,
                                 save_path=OUT_FIXED_HAP) -> pd.DataFrame:
        id_hap = pd.read_csv(id_hap_path, sep='\t')
        validate_input(id_hap, ["ID", "Haplogroup"])

        corr_map = pd.read_csv(correction_path).set_index("Origin")["Revised"]
        id_hap["Haplogroup"] = id_hap["Haplogroup"].replace(corr_map.to_dict())

        id_hap.to_csv(save_path, sep='\t', index=False)
        return id_hap

    # 第二步：构建 PhyloTree 追溯链
    def build_parent_map(tree_path=FILE_PHYLOTREE) -> dict[str, str | None]:
        parent_map: dict[str, str | None] = {}
        stack: list[tuple[int, str]] = []
        with open(tree_path, encoding='utf-8') as f:
            for raw in f:
                if not raw.strip():
                    continue
                level = raw.count('\t')
                hap = raw.strip()
                while stack and stack[-1][0] >= level:
                    stack.pop()
                parent = stack[-1][1] if stack else None
                parent_map[hap] = parent
                stack.append((level, hap))
        return parent_map

    def build_lineage_functions(parent_map):
        @lru_cache(maxsize=None)
        def lineage(hap: str) -> tuple[str, ...]:
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
        id_hap["lineage"] = id_hap["Haplogroup"].map(lineage_func)
        max_len = id_hap["lineage"].str.len().max()

        level_cols_fwd = [f"Level_{i}" for i in range(max_len)]
        level_cols_rev = [f"Level_{i}" for i in range(max_len-1, -1, -1)]

        def pad(seq: tuple[str, ...]) -> list[str]:
            seq = list(seq)
            return seq + [""] * (max_len - len(seq))

        id_hap[level_cols_fwd] = id_hap["lineage"].apply(lambda x: pd.Series(pad(x)))
        id_hap[level_cols_rev] = id_hap["lineage"].apply(lambda x: pd.Series(pad(tuple(reversed(x)))))

        id_hap[["ID"] + level_cols_fwd].to_csv(forward_file, sep='\t', index=False)
        id_hap[["ID"] + level_cols_rev].to_csv(reverse_file, sep='\t', index=False)

        missed = id_hap[id_hap["lineage"].str.len() == 0][["ID", "Haplogroup"]]
        if not missed.empty:
            missed.to_csv(not_found_file, sep='\t', index=False)

    # 第三步：匹配目标列表
    def match_to_target(reverse_file=OUT_REVERSED_LEVEL,
                        target_path=FILE_TARGET,
                        rough_fix_path=FILE_ROUGH_FIX,
                        out_final=OUT_FINAL,
                        out_unmatched=OUT_UNMATCHED) -> None:
        df_rev = pd.read_csv(reverse_file, sep='\t')
        validate_input(df_rev, ["ID"])

        target_set = set(pd.read_csv(target_path, header=None, names=["Hap"])['Hap'])
        rough_map = pd.read_csv(rough_fix_path, sep='\t', header=None,
                                names=['Original', 'Correction']).set_index("Original")["Correction"]

        level_cols = df_rev.columns[1:]

        matched, unmatched = [], []
        for _, row in df_rev.iterrows():
            seq = [h for h in row[level_cols] if pd.notna(h) and str(h).strip()]
            hit = next((h for h in seq if h in target_set), None)
            if hit:
                matched.append((row["ID"], rough_map.get(hit, hit)))
            else:
                fall_back = seq[0] if seq else ""
                unmatched.append((row["ID"], fall_back))

        pd.DataFrame(matched, columns=["ID", "Haplogroup_PCA"]).to_csv(out_final, sep='\t', index=False)
        pd.DataFrame(unmatched, columns=["ID", "Haplogroup"]).to_csv(out_unmatched, sep='\t', index=False)

        if unmatched and all(h.startswith("L") for _, h in unmatched if h):
            df_un = pd.read_csv(out_unmatched, sep='\t')
            df_un["Haplogroup_PCA"] = df_un["Haplogroup"].str[:2]
            pd.concat([pd.read_csv(out_final, sep='\t'), df_un[["ID", "Haplogroup_PCA"]]]
                      ).to_csv(out_final, sep='\t', index=False)

    # ==================== 执行流程 ====================
    try:
        df_hap = apply_primary_correction()
        lineage = build_lineage_functions(build_parent_map())
        output_level_tables(df_hap, lineage)
        match_to_target()
        print("✅ YuChunLi流程完毕")
    except Exception as e:
        print(f"💥 YuChunLi流程中断: {e}", file=sys.stderr)
        raise


# ============================================================================
# 步骤9-11：LLT流程（来自笔记本单元格9 - 规范单倍群）
# ============================================================================
def run_llt_workflow(id_hap_path: str, conf_dir: str, temp_dir: str):
    """
    LLT方案流程
    代码直接来自笔记本单元格9，只改路径部分
    """
    print("\n" + "="*60)
    print("🔴 LLT 流程")
    print("="*60)
    
    # 配置路径（只改这部分）
    OUTPUT_DIR = Path(temp_dir)
    SRC_DIR = Path(conf_dir)
    
    FILE_ID_HAP = Path(id_hap_path)
    FILE_CORRECTION = SRC_DIR / "Haplogrep单倍群分型订正表.csv"
    FILE_PHYLOTREE = SRC_DIR / "线粒体单倍群phylotree(version17)2025年3月12日.txt"
    FILE_TARGET = SRC_DIR / "目标_LLT.txt"
    FILE_ROUGH_FIX = SRC_DIR / "错误纠正_LLT.txt"
    
    OUT_FIXED_HAP = OUTPUT_DIR / "订正之后的单倍群名称_LLT.txt"
    OUT_FORWARD_LEVEL = OUTPUT_DIR / "正序等级_LLT.txt"
    OUT_REVERSED_LEVEL = OUTPUT_DIR / "逆序等级_LLT.txt"
    OUT_NOT_FOUND = OUTPUT_DIR / "没有查询到请核实_LLT.txt"
    OUT_FINAL = OUTPUT_DIR / "最终_LLT.txt"
    OUT_UNMATCHED = OUTPUT_DIR / "无法处理_LLT.txt"
    
    # ==================== 以下代码直接来自笔记本，无任何修改 ====================
    def validate_input(df: pd.DataFrame, required: list[str]) -> None:
        if not set(required) <= set(df.columns):
            raise ValueError(f"输入缺少必要列: {required}")

    # 第一步：纠正 haplogroup 名称
    def apply_primary_correction(id_hap_path=FILE_ID_HAP,
                                 correction_path=FILE_CORRECTION,
                                 save_path=OUT_FIXED_HAP) -> pd.DataFrame:
        id_hap = pd.read_csv(id_hap_path, sep='\t')
        validate_input(id_hap, ["ID", "Haplogroup"])

        corr_map = pd.read_csv(correction_path).set_index("Origin")["Revised"]
        id_hap["Haplogroup"] = id_hap["Haplogroup"].replace(corr_map.to_dict())

        id_hap.to_csv(save_path, sep='\t', index=False)
        return id_hap

    # 第二步：构建 PhyloTree 追溯链
    def build_parent_map(tree_path=FILE_PHYLOTREE) -> dict[str, str | None]:
        parent_map: dict[str, str | None] = {}
        stack: list[tuple[int, str]] = []
        with open(tree_path, encoding='utf-8') as f:
            for raw in f:
                if not raw.strip():
                    continue
                level = raw.count('\t')
                hap = raw.strip()
                while stack and stack[-1][0] >= level:
                    stack.pop()
                parent = stack[-1][1] if stack else None
                parent_map[hap] = parent
                stack.append((level, hap))
        return parent_map

    def build_lineage_functions(parent_map):
        @lru_cache(maxsize=None)
        def lineage(hap: str) -> tuple[str, ...]:
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
        id_hap["lineage"] = id_hap["Haplogroup"].map(lineage_func)
        max_len = id_hap["lineage"].str.len().max()

        level_cols_fwd = [f"Level_{i}" for i in range(max_len)]
        level_cols_rev = [f"Level_{i}" for i in range(max_len-1, -1, -1)]

        def pad(seq: tuple[str, ...]) -> list[str]:
            seq = list(seq)
            return seq + [""] * (max_len - len(seq))

        id_hap[level_cols_fwd] = id_hap["lineage"].apply(lambda x: pd.Series(pad(x)))
        id_hap[level_cols_rev] = id_hap["lineage"].apply(lambda x: pd.Series(pad(tuple(reversed(x)))))

        id_hap[["ID"] + level_cols_fwd].to_csv(forward_file, sep='\t', index=False)
        id_hap[["ID"] + level_cols_rev].to_csv(reverse_file, sep='\t', index=False)

        missed = id_hap[id_hap["lineage"].str.len() == 0][["ID", "Haplogroup"]]
        if not missed.empty:
            missed.to_csv(not_found_file, sep='\t', index=False)

    # 第三步：匹配目标列表
    def match_to_target(reverse_file=OUT_REVERSED_LEVEL,
                        target_path=FILE_TARGET,
                        rough_fix_path=FILE_ROUGH_FIX,
                        out_final=OUT_FINAL,
                        out_unmatched=OUT_UNMATCHED) -> None:
        df_rev = pd.read_csv(reverse_file, sep='\t')
        validate_input(df_rev, ["ID"])

        target_set = set(pd.read_csv(target_path, header=None, names=["Hap"])['Hap'])
        rough_map = pd.read_csv(rough_fix_path, sep='\t', header=None,
                                names=['Original', 'Correction']).set_index("Original")["Correction"]

        level_cols = df_rev.columns[1:]

        matched, unmatched = [], []
        for _, row in df_rev.iterrows():
            seq = [h for h in row[level_cols] if pd.notna(h) and str(h).strip()]
            hit = next((h for h in seq if h in target_set), None)
            if hit:
                matched.append((row["ID"], rough_map.get(hit, hit)))
            else:
                fall_back = seq[0] if seq else ""
                unmatched.append((row["ID"], fall_back))

        pd.DataFrame(matched, columns=["ID", "Haplogroup_PCA"]).to_csv(out_final, sep='\t', index=False)
        pd.DataFrame(unmatched, columns=["ID", "Haplogroup"]).to_csv(out_unmatched, sep='\t', index=False)

        if unmatched and all(h.startswith("L") for _, h in unmatched if h):
            df_un = pd.read_csv(out_unmatched, sep='\t')
            df_un["Haplogroup_PCA"] = df_un["Haplogroup"].str[:2]
            pd.concat([pd.read_csv(out_final, sep='\t'), df_un[["ID", "Haplogroup_PCA"]]]
                      ).to_csv(out_final, sep='\t', index=False)

    # ==================== 执行流程 ====================
    try:
        df_hap = apply_primary_correction()
        lineage = build_lineage_functions(build_parent_map())
        output_level_tables(df_hap, lineage)
        match_to_target()
        print("✅ LLT流程完毕")
    except Exception as e:
        print(f"💥 LLT流程中断: {e}", file=sys.stderr)
        raise


# ============================================================================
# 步骤12：合并结果（来自笔记本单元格11）
# ============================================================================
def merge_results(temp_dir: str, data_dir: str):
    """
    合并YuChunLi和LLT的结果
    """
    print("\n" + "="*60)
    print("🟢 合并结果")
    print("="*60)
    
    try:
        temp_path = Path(temp_dir)
        data_path = Path(data_dir)
        data_path.mkdir(parents=True, exist_ok=True)
        
        df_ID_hap = pd.read_csv(temp_path / '订正之后的单倍群名称_YuChunLi.txt', sep='\t')
        df_LLT = pd.read_csv(temp_path / '最终_LLT.txt', sep='\t').rename(columns={'Haplogroup_PCA': 'Haplogroup_LLT'})
        df_YCL = pd.read_csv(temp_path / '最终_YuChunLi.txt', sep='\t').rename(columns={'Haplogroup_PCA': 'Haplogroup_YuChunLi'})
        dF_最终 = df_ID_hap.merge(df_LLT, on='ID', how='left')
        dF_最终 = dF_最终.merge(df_YCL, on='ID', how='left')
        dF_最终.to_csv(data_path / '全部单倍群整理.csv', index=False)
        
        print(f"✅ 合并完成：{data_path / '全部单倍群整理.csv'}")
    except Exception as e:
        print(f"❌ 合并失败: {e}", file=sys.stderr)
        raise


def main():
    """命令行入口"""
    parser = argparse.ArgumentParser(
        description="单倍群整理完整流程：按照Jupyter笔记本执行"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入Excel文件"
    )
    parser.add_argument(
        "--config-dir", "-c",
        required=True,
        help="配置文件目录（包含纠正表、PhyloTree等）"
    )
    parser.add_argument(
        "--temp-dir", "-t",
        required=True,
        help="临时文件目录"
    )
    parser.add_argument(
        "--data-dir", "-d",
        required=True,
        help="数据输出目录"
    )
    
    args = parser.parse_args()
    
    try:
        print("\n" + "="*60)
        print("🚀 单倍群整理完整流程")
        print("="*60)
        
        # 步骤1-5：提取原始单倍群
        print("\n📊 [步骤1-5] 提取原始单倍群...")
        temp_path = Path(args.temp_dir)
        temp_path.mkdir(parents=True, exist_ok=True)
        
        df_raw = extract_raw_haplogroups(args.input)
        
        id_hap_path = temp_path / "ID_Hap.txt"
        df_raw.to_csv(id_hap_path, sep='\t', index=False)
        print(f"✅ 提取完成：{len(df_raw)} 条记录")
        print(f"   保存到：{id_hap_path}")
        
        # 步骤6-8：YuChunLi流程
        run_yuchunli_workflow(str(id_hap_path), args.config_dir, args.temp_dir)
        
        # 步骤9-11：LLT流程
        run_llt_workflow(str(id_hap_path), args.config_dir, args.temp_dir)
        
        # 步骤12：合并结果
        merge_results(args.temp_dir, args.data_dir)
        
        print("\n" + "="*60)
        print("✅ 全部流程运行完毕！")
        print("="*60 + "\n")
        
        sys.exit(0)
    except Exception as e:
        print(f"❌ 错误：{e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
