#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/12/15 15:25
# @Last Modified by:   Ming
# @Last Modified time: 2022/12/15 15:25
import logging
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--pairedkmers',
              required=True,
              type=click.Path(),
              help="The pairdekmers dir of the pipe")
@click.option('--chromosomes',
              required=True,
              type=click.Path(),
              help="The chromosome name list sep by ','")
@click.option('--name',
              required=True,
              type=click.Path(),
              help="The Latin name")
@click.option('--tolerance',
              default=80,
              show_default=True,
              type=float,
              help="包容性阈值")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The out put bed file name")
def main(pairedkmers, chromosomes, name, tolerance, out):
    """
    合并并处理paired-kmers生成的结果文件
    """
    d_pairedkmers = Path(pairedkmers).absolute()
    l_chromosome = chromosomes.strip().split(',')
    logger.info(f"Start to Deal the paired-kmers result")
    f_out = Path(out).absolute()
    flag = 1
    with open(f_out, 'w') as OUT:
        for chromosome in l_chromosome:
            f_in = d_pairedkmers.joinpath(f"{chromosome}/kmer.bed")
            with open(f_in, 'r') as IN:
                for line in IN:
                    if line.startswith('#'):
                        pass
                    else:
                        arr = line.strip().split("\t")
                        if arr[0] == chromosome and float(arr[4]) >= tolerance:
                            arr[1] = int(arr[1]) - len(arr[5])
                            arr[2] = int(arr[2]) + len(arr[6])
                            print(*[arr[0], arr[1], arr[2], f"{name}-pkmers{flag}"], sep="\t", file=OUT)
                            flag += 1


if __name__ == "__main__":
    main()
