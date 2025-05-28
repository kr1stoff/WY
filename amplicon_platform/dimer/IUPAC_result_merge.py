import click
import re


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input', required=True, type=click.Path(), help="输入")
@click.option('-o', '--output', default='./out.txt', type=click.Path(), help="输出")
def main(input, output):
    # 打开文件并逐行读取数据
    with open(input, 'r') as file:
        data_dict = {}  # 创建一个字典用于存储数据
        for line in file:
            line = line.strip()
            arr = line.split('\t')
            dimer = arr[1]
            score = arr[2]
            # S0110-10__F__in[1]::S0110-11__R__in[2]
            pattern = r'(\S+)\[\d+\]::(\S+)\[\d+\]'
            match = re.match(pattern, dimer)
            if match:
                p1, p2 = match.groups()
                primer_group = f'{p1}::{p2}'
                if primer_group not in data_dict or score > data_dict[primer_group][0]:
                    data_dict[primer_group] = (score, line)  # 保留最大的C值

    # 输出结果
    with open(output, 'w') as output_file:
        for key, value in data_dict.items():
            arr = value[1].split('\t')
            arr[1] = key
            out = '\t'.join(arr)
            print(out, file=output_file)


if __name__ == '__main__':
    main()
