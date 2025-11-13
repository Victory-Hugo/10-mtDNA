#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
功能：根据样本原始单倍群的首字母，重新调整等级表
输入：ID_Hap.txt 和 正序等级-补齐至下游.txt
输出：正序等级-补齐至下游-按首字母.txt

逻辑：
对于每个样本，获取其原始单倍群（如 M7b, B4a1c2 等）
提取首字母（如 M, B）
在该样本的正序等级-补齐至下游.txt 中的一行中，找到首个以该首字母开头的单倍群
将该单倍群及其之后的所有单倍群作为新的 Level_0, Level_1, ...
"""

import pandas as pd
from pathlib import Path
import sys

def extract_first_letter(haplogroup: str) -> str:
    """提取单倍群的首字母"""
    if not haplogroup or pd.isna(haplogroup):
        return ""
    return str(haplogroup)[0]

def find_start_level(row_values, first_letter: str) -> int:
    """
    在一行数据中找到首个以指定首字母开头的单倍群的索引
    特殊处理：如果首字母是 L，则从 MT-MRCA 开始（索引0）
    返回该索引，如果找不到则返回 -1
    """
    # 如果首字母是 L，统一从 MT-MRCA 开始（即不删除任何列）
    if first_letter == 'L':
        return 0
    
    for i, val in enumerate(row_values):
        if pd.notna(val) and str(val).strip():
            if str(val)[0] == first_letter:
                return i
    return -1

def process_by_first_letter(id_hap_file: Path, 
                           forward_level_file: Path, 
                           output_file: Path) -> None:
    """
    处理数据，按首字母重新调整等级
    """
    # 读取原始单倍群信息
    id_hap_df = pd.read_csv(id_hap_file, sep='\t')
    id_hap_dict = dict(zip(id_hap_df['ID'], id_hap_df['Haplogroup']))
    
    # 读取正序等级-补齐至下游文件
    level_df = pd.read_csv(forward_level_file, sep='\t')
    
    # 创建新的结果 DataFrame
    result_rows = []
    
    for idx, row in level_df.iterrows():
        sample_id = row['ID']
        
        # 获取样本的原始单倍群
        original_hap = id_hap_dict.get(sample_id, "")
        first_letter = extract_first_letter(original_hap)
        
        if not first_letter:
            # 如果找不到原始单倍群，保持原样
            result_rows.append(row.tolist())
            continue
        
        # 获取该样本的所有等级值
        level_values = row.iloc[1:].tolist()  # 排除 ID 列
        
        # 找到首个以该首字母开头的单倍群的位置
        start_idx = find_start_level(level_values, first_letter)
        
        if start_idx == -1:
            # 如果找不到以该首字母开头的单倍群，保持原样
            result_rows.append(row.tolist())
            continue
        
        # 从 start_idx 开始截取数据
        new_values = level_values[start_idx:]
        
        # 如果数据太短，填充空值使其长度与原来一致
        while len(new_values) < len(level_values):
            new_values.append(new_values[-1] if new_values else "")
        
        # 重新组合：ID + 截取后的值
        new_row = [sample_id] + new_values
        result_rows.append(new_row)
    
    # 创建新 DataFrame
    column_names = level_df.columns.tolist()
    result_df = pd.DataFrame(result_rows, columns=column_names)
    
    # 保存结果
    output_file.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"✅ 完成: {output_file}")
    print(f"   - ID信息文件: {id_hap_file}")
    print(f"   - 正序等级文件: {forward_level_file}")
    print(f"   - 输出文件: {output_file}")
    print(f"   - 处理的样本数: {len(result_df)}")

if __name__ == '__main__':
    id_hap_file = Path("/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/2-母系专用/example/ID_Hap.txt")
    forward_level_file = Path("/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/2-母系专用/example/output/正序等级-补齐至下游.txt")
    output_file = Path("/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/2-母系专用/example/output/正序等级-补齐至下游-按首字母.txt")
    
    try:
        process_by_first_letter(id_hap_file, forward_level_file, output_file)
    except Exception as e:
        print(f"❌ 错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
