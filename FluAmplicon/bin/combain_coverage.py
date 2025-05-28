#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/24 12:37
# @Last Modified by:   Ming
# @Last Modified time: 2022/6/24 12:37
import logging
import os.path

import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-d', '--dirin',
              required=True,
              type=click.Path(),
              help="The dir contain the coverage info")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put file")
def main(dirin, out):
    """
    合并深度统计结果
    """
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    res = []
    for segment in segments:
        f_in = os.path.join(dirin, f'{segment}.bam_stats.txt')
        if os.path.exists(f_in):
            with open(f_in, 'r') as IN:
                next(IN)
                for line in IN:
                    arr = line.strip().split('\t')
                    res.append(arr)
        else:
            res.append(["NA"] * 6)

    with open(out, 'w') as OUT:
        header = ["片段", "平均深度", "覆盖度", "深度≥10x", "深度≥30x", "深度≥100x", "均一性"]
        print(*header, sep="\t", file=OUT)
        for i, j in zip(segments, res):
            arr = [i] + j
            print(*arr, sep='\t', file=OUT)


if __name__ == "__main__":
    main()
