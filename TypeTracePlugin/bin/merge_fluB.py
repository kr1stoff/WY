#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/7/12 10:48
# @Last Modified by:   Ming
# @Last Modified time: 2022/7/12 10:48
import logging

import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('--flu',
              required=True,
              type=click.Path(),
              help="The result get by LABEL irma-FLU-v2")
@click.option('--segments',
              required=True,
              help="The result get by LABEL segment module(sep by ,)")
@click.option('--out',
              required=True,
              help="The final merge result")
def main(flu, segments, out):
    """
    合并LABEL的乙流分型结果中irma-FLU-v2模块和其它各个segment模块的结果
    """
    logging.info(f"Read the file {flu}")
    res = {}
    with open(flu, 'r') as IN:
        next(IN)
        for line in IN:
            arr = [i.strip() for i in line.strip().split('\t')]
            res.setdefault(arr[0], ["-", "-"])
            if arr[1] != "UNRECOGNIZABLE":
                res[arr[0]][0] = arr[1]

    for segment in segments.strip().split(','):
        logging.info(f"Parse the file {segment}")
        with open(segment, 'r') as IN:
            next(IN)
            for line in IN:
                arr = [i.strip() for i in line.strip().split('\t')]
                if arr[1] != "UNRECOGNIZABLE":
                    res[arr[0]][1] = arr[1]

    with open(out, 'w') as OUT:
        print(*["序列名称", "片段预测", "分型预测"], sep="\t", file=OUT)
        for k, v in res.items():
            print(*[k, *v], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
