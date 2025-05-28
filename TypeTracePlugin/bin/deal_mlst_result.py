#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/7/5 14:41
# @Last Modified by:   Ming
# @Last Modified time: 2022/7/5 14:41
import logging

import click


#### Some Function
def get_min(content):
    """
    获取最小的等位基因编号
    """
    if '-' in content:
        return '-'
    elif ',' in content or '?' in content:
        res = str(min([int(i.replace('?', '')) for i in content.strip().split(',')]))
        return res
    else:
        return content


##################


__version__ = '1.1.3'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('--input',
              required=True,
              type=click.Path(),
              help="The input mlst result")
@click.option('--scheme',
              required=True,
              type=click.Path(),
              help="The species scheme file")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The output file")
def main(input, scheme, out):
    """
    处理mlst的结果

    目的：
    1. 部分结果文件的头和内容列数不一样，将缺失的内容补为"-"
    2. 去除基因分型序号的中的？
    3. 有多个等位基因的，取最小值
    """
    logging.info(f"Read the ori mlst result")
    with open(input, 'r') as IN:
        header = next(IN).strip().split("\t")
        arr = next(IN).strip().split('\t')
        if len(header) > len(arr):
            arr += ["-" for i in range(len(header) - len(arr))]
        elif len(header) < len(arr):
            arr = arr[:len(header)]
        else:
            pass

        new_header = header[2:]
        n_arr = arr[2:]
        if n_arr[0] == "-":
            logging.info(f"Get scheme info")
            num_gene = len(new_header) - 1
            info = {}
            with open(scheme, 'r') as IN:
                next(IN)
                for line in IN:
                    arr = line.strip().split("\t")
                    key = tuple(arr[1:num_gene + 1])
                    info[key] = arr[0]
            tmp = [get_min(i) for i in n_arr[1:]]
            if tuple(tmp) in info:
                n_arr[0] = info[tuple(tmp)]
        else:
            pass

    logging.info(f"Out put the result to {out}")
    with open(out, 'w') as OUT:
        print(*new_header, sep="\t", file=OUT)
        print(*n_arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
