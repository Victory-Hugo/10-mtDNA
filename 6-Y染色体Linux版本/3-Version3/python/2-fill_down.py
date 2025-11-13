#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
功能：将正序等级文件中每一行最右侧的单倍群填充至最后一列，消除空值
输入：正序等级.txt
输出：正序等级-补齐至下游.txt
"""

import pandas as pd
from pathlib import Path
import sys
import argparse

def fill_forward(input_file: Path, output_file: Path) -> None:
    """
    读取输入文件，将每行最右侧的非空值向右填充到最后一列
    """
    # 读取文件
    df = pd.read_csv(input_file, sep='\t')
    
    # 保存 ID 列
    id_col = df.iloc[:, 0]
    
    # 对数据部分（除了 ID 列）进行前向填充
    data_cols = df.iloc[:, 1:]
    
    # 对每一行进行前向填充（ffill），将右侧的非空值填充至末尾的空值
    # 由于我们需要向右填充，使用 bfill 在轴向1（列向）上
    # 但 bfill 是从右到左填充，所以需要反转后再处理
    data_cols_filled = data_cols.ffill(axis=1)
    
    # 重新组合 ID 列和填充后的数据列
    result = pd.concat([id_col, data_cols_filled], axis=1)
    
    # 保存结果
    output_file.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_file, sep='\t', index=False)
    
    print(f"✅ 完成: {output_file}")
    print(f"   - 输入文件: {input_file}")
    print(f"   - 输出文件: {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='将正序等级文件中每一行最右侧的单倍群填充至最后一列',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input-file', '-i', required=True,
                        help='输入文件路径（正序等级.txt）')
    parser.add_argument('--output-file', '-o', required=True,
                        help='输出文件路径（正序等级-补齐至下游.txt）')
    
    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_file = Path(args.output_file)
    
    try:
        fill_forward(input_file, output_file)
    except Exception as e:
        print(f"❌ 错误: {e}", file=sys.stderr)
        sys.exit(1)
