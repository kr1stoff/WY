#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/19 16:00
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/19 16:00
import logging
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-s', '--specificity',
              required=True,
              type=click.Path(),
              help="The specificity merge result file")
@click.option('-i', '--inclusive',
              required=True,
              type=click.Path(),
              help="The inclusive merge result file")
@click.option('--info',
              required=True,
              type=click.Path(),
              help="The pPCR pipe result or after get the degenerate")
@click.option('--species',
              required=True,
              type=click.Path(),
              help="The species id")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put result file")
def main(specificity, inclusive, info, species, out):
    """
    合并特异性以及包容性评估的结果，生成最终的文件
    """
    f_spe = Path(specificity).absolute()
    f_inc = Path(inclusive).absolute()
    f_info = Path(info).absolute()
    f_out = Path(out).absolute()

    logger.info(f"Get the primer info from {f_info}")
    info_primer = {}
    with open(f_info, 'r') as IN:
        header_primer = next(IN).strip().split("\t")
        for line in IN:
            arr = line.strip().split("\t")
            info_primer[arr[0]] = arr

    logger.info(f"Get the specificity info from {f_spe}")
    info_spe = {}
    with open(f_spe, 'r') as IN:
        header_spe = [f"{i}_Specificity" for i in next(IN).strip().split("\t")[1:]]
        for line in IN:
            arr = line.strip().split("\t")
            info_spe[arr[0]] = arr[1:]

    logger.info(f"Get the inclusiver info from {f_inc}")
    info_inc = {}
    with open(f_inc, 'r') as IN:
        header_inc = [f"{i}_Inclusive" for i in next(IN).strip().split("\t")[1:]]
        for line in IN:
            arr = line.strip().split("\t")
            info_inc[arr[0]] = arr[1:]

    logger.info(f"Out put the result to file {f_out}")
    with open(f_out, 'w') as OUT:
        header = ["Species"] + header_primer + header_spe + header_inc
        print(*header, sep="\t", file=OUT)
        for i, j in info_primer.items():
            arr = [species, *j, *info_spe[i], *info_inc[i]]
            print(*arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
