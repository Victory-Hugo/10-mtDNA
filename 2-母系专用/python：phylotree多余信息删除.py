# 这个脚本的作用是将excel全部线粒体单倍群发育树转化成可以被直接使用的txt文件
# 定义文件路径
file_path = 'C:/Users/victo/Desktop/新建 Text Document.txt'

# 读取文件并逐行处理
processed_lines = []
with open(file_path, 'r', encoding='utf-8') as file:
    for line in file:
        # 去掉行首的制表符，以识别第一个非空字符
        stripped_line = line.lstrip('\t')
        
        # 检查行是否包含任何非空字符
        if stripped_line:
            # 找到第一个非空字符的位置
            non_whitespace_index = line.find(stripped_line[0])
            # 分割去掉空格后的行并获取第一个单词
            first_word = stripped_line.split()[0] if stripped_line.split() else ""
            first_word_end = non_whitespace_index + len(first_word)
            
            # 在剩余部分查找第一个制表符的位置
            tab_index = line.find('\t', first_word_end)
            
            # 如果找到了制表符，去掉制表符之后的内容
            if tab_index != -1:
                processed_line = line[:tab_index].rstrip()
            else:
                processed_line = line.rstrip()

            processed_lines.append(processed_line)

# 定义输出文件的路径
output_file_path = 'C:/Users/victo/Desktop/processed_text.txt'

# 将处理后的行写入输出文件
with open(output_file_path, 'w', encoding='utf-8') as output_file:
    output_file.write('\n'.join(processed_lines))

# 输出文件路径
output_file_path
