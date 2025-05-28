import ruamel.yaml

# 创建一个有序字典，按照指定的顺序排列键值对
data = {
    'nt_database': 1,
    'ptype': 'Bacteria',
    'resion_type': 'gene'
}

# 创建一个列表，表示键值对的顺序
key_order = ['nt_database', 'ptype', 'resion_type']

# 创建注释字典
comments = {
    'nt_database': 'taxid',
}

# 构建 YAML 文件内容
yaml_content = ""
for key in key_order:
    if key in comments:
        yaml_content += f"# {comments[key]}\n"
    yaml_content += f"{key}: {data[key]}\n"

# 将内容写入 YAML 文件
with open('output.yaml', 'w') as yaml_file:
    yaml_file.write(yaml_content)
