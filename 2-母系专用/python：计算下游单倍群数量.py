# 文件路径
haplogroups_file = 'C:/Users/victo/Desktop/新建 Text Document.txt'
phylogroup_file = 'F:/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/线粒体单倍群phylotree(version17).txt'
output_file = 'C:/Users/victo/Desktop/haplogroup_counts.txt'

def read_haplogroups(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    return lines

def find_downstream_haplogroups(haplogroup, lines):
    downstream_haplogroups = []
    haplogroup_level = None
    capture = False

    for line in lines:
        stripped_line = line.strip()
        if not stripped_line:
            continue
        
        current_level = len(line) - len(line.lstrip('\t'))
        if haplogroup_level is None and stripped_line.startswith(haplogroup):
            haplogroup_level = current_level
            capture = True
            continue
        
        if capture:
            if current_level > haplogroup_level:
                downstream_haplogroups.append(stripped_line)
            elif current_level == haplogroup_level:
                break

    return downstream_haplogroups

def count_haplogroups(haplogroups_list, target_haplogroup, phylogroup_lines):
    count = 0
    visited = set()

    def recursive_count(haplogroup):
        nonlocal count
        if haplogroup in visited:
            return
        visited.add(haplogroup)
        haplogroup_count = haplogroups_list.count(haplogroup)
        count += haplogroup_count
        downstream_haplogroups = find_downstream_haplogroups(haplogroup, phylogroup_lines)
        for downstream in downstream_haplogroups:
            recursive_count(downstream)
    
    recursive_count(target_haplogroup)
    return count

def main():
    # 读取目标单倍群
    haplogroups = read_haplogroups(haplogroups_file)
    haplogroups_list = [haplo.strip() for haplo in haplogroups]
    
    # 读取所有单倍群数据
    phylogroup_lines = read_haplogroups(phylogroup_file)
    
    # 查找每个单倍群的下游单倍群并统计数量
    results = {}
    for haplogroup in haplogroups_list:
        if haplogroup not in results:
            count = count_haplogroups(haplogroups_list, haplogroup, phylogroup_lines)
            downstream = find_downstream_haplogroups(haplogroup, phylogroup_lines)
            results[haplogroup] = {
                'count': count,
                'downstream': downstream
            }
    
    # 输出结果到文件
    with open(output_file, 'w', encoding='utf-8') as f:
        for haplogroup, data in results.items():
            f.write(f'{haplogroup}\t{data["count"]}\n')

    print(f'结果已保存到 {output_file}')

if __name__ == '__main__':
    main()
