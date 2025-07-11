import csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def parse_haplogroup_file(lines):
    haplogroups = []
    for line in lines:
        stripped_line = line.strip()
        if stripped_line:
            level = line.count('\t')
            haplogroups.append((level, stripped_line))
    return haplogroups

def find_correct_upstream_haplogroups(haplogroup_name, haplogroups):
    upstream_haplogroups = []
    current_haplogroup = None
    for level, haplogroup in reversed(haplogroups):
        if haplogroup == haplogroup_name:
            current_haplogroup = haplogroup
            current_level = level
            upstream_haplogroups.append((level, haplogroup))
        elif current_haplogroup and level < current_level:
            upstream_haplogroups.insert(0, (level, haplogroup))
            current_level = level
    return upstream_haplogroups

def process_haplogroup(haplogroup_name, haplogroups):
    correct_upstream_haplogroups = find_correct_upstream_haplogroups(haplogroup_name, haplogroups)
    if not correct_upstream_haplogroups:
        return haplogroup_name, None, None
    haplogroup_list = [hap for _, hap in correct_upstream_haplogroups]
    result_line = "\t".join(haplogroup_list)
    reversed_result_line = "\t".join(haplogroup_list[::-1])
    return haplogroup_name, result_line, reversed_result_line

def process_wrapper(args):
    haplogroup_name, haplogroups = args
    return process_haplogroup(haplogroup_name, haplogroups)

def process_haplogroups(
    input_file_path: Path,
    output_file_path: Path,
    reversed_output_file_path: Path,
    not_found_file_path: Path,
    haplogroup_file_path: Path
):
    # 读取单倍群文件
    with haplogroup_file_path.open('r', encoding='utf-8') as file:
        refined_lines = file.readlines()
    haplogroups = parse_haplogroup_file(refined_lines)
    
    # 读取输入文件，确保包含 ID 和 Haplogroup 列
    with input_file_path.open('r', encoding='utf-8') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        fieldnames = reader.fieldnames
        if not fieldnames or 'ID' not in fieldnames or 'Haplogroup' not in fieldnames:
            raise ValueError("输入文件必须包含 'ID' 和 'Haplogroup' 两列。")
        rows = list(reader)
    
    # 去重单倍群名称，避免重复计算
    unique_haplogroup_names = list({row['Haplogroup'] for row in rows})
    args = [(name, haplogroups) for name in unique_haplogroup_names]
    
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_wrapper, args))
    
    results_dict = {name: (result_line, reversed_result_line) for name, result_line, reversed_result_line in results}
    not_found_entries = []
    
    # 计算最大层级数（用于输出对齐）
    valid_results = [result_line for result_line, _ in results_dict.values() if result_line is not None]
    if valid_results:
        max_level = max(len(result_line.split("\t")) for result_line in valid_results)
    else:
        max_level = 0
    
    with output_file_path.open('w', encoding='utf-8') as output_file, \
         reversed_output_file_path.open('w', encoding='utf-8') as reversed_output_file:
        header = ["ID"] + [f"Level_{i}" for i in range(max_level)]
        output_file.write("\t".join(header) + "\n")
        reversed_header = ["ID"] + [f"Level_{i}" for i in range(max_level-1, -1, -1)]
        reversed_output_file.write("\t".join(reversed_header) + "\n")
    
        for row in rows:
            id_value = row['ID']
            haplogroup_name = row['Haplogroup']
            if haplogroup_name in results_dict:
                result_line, reversed_result_line = results_dict[haplogroup_name]
                if result_line:
                    levels = result_line.split("\t")
                    reversed_levels = reversed_result_line.split("\t")
                    if len(levels) < max_level:
                        levels += [""] * (max_level - len(levels))
                    if len(reversed_levels) < max_level:
                        reversed_levels += [""] * (max_level - len(reversed_levels))
                    output_file.write("\t".join([id_value] + levels) + "\n")
                    reversed_output_file.write("\t".join([id_value] + reversed_levels) + "\n")
                else:
                    not_found_entries.append((id_value, haplogroup_name))
            else:
                not_found_entries.append((id_value, haplogroup_name))
    
    if not_found_entries:
        with not_found_file_path.open('w', encoding='utf-8') as not_found_file:
            not_found_file.write("ID\tHaplogroup\n")
            for id_value, haplogroup_name in set(not_found_entries):
                not_found_file.write(f"{id_value}\t{haplogroup_name}\n")

if __name__ == "__main__":
    # 示例：所有变量均传入绝对路径
    process_haplogroups(
        input_file_path=Path("需要查询的文件.txt"),
        output_file_path=Path("正序等级.txt"),
        reversed_output_file_path=Path("逆序等级.txt"),
        not_found_file_path=Path("没有查询到请核实.txt"),
        haplogroup_file_path=Path("线粒体单倍群phylotree(version17)2025年3月12日.txt")
    )