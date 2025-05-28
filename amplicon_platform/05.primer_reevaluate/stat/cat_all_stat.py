import os
import click


def get_path(id, file_path, file_type):
    file_type_new = file_type.replace('$id', id)
    input_file_path = os.path.join(file_path, file_type_new)
    return input_file_path


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_file', required=True, type=click.Path(), help="表格")
@click.option('-f', '--file_path', required=True, type=click.Path(), help="文件路径")
@click.option('-t', '--file_type', default='$id/$id.all.stat', type=str, help="类型")
@click.option('-o', '--out', required=True, type=click.Path(), help="输出")
def main(input_file, file_path, file_type, out):
    # 打开输出文件以准备写入数据
    with open(out, "w") as outfile:
        # 打开输入文件以逐行读取
        with open(input_file, "r") as IN:
            for line in IN:
                id = line.strip().split('\t')[0]
                input_file_path = get_path(id, file_path, file_type)
                if os.path.exists(input_file_path) and os.path.getsize(input_file_path) != 0:
                    # 打开输入文件以逐行读取并筛选
                    with open(input_file_path, "r") as inputfile:
                        for input_line in inputfile:
                            if not input_line.startswith("#"):
                                outfile.write(input_line)


if __name__ == '__main__':
    main()
