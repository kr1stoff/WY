#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/1/28 15:07
# @Last Modified by:   Ming
# @Last Modified time: 2022/1/28 15:07
import logging

import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-f', '--fname',
              required=True,
              type=click.Path(),
              help="The DADA2 output stat file name")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put result for plot")
def main(fname, out):
    """
    Get the data for plot for denois stat
    """
    info = {}
    logging.info(f"Read the input file")
    with open(fname, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split('\t')
            info[arr[0]] = [int(i) for i in [arr[1], arr[2], arr[4], arr[5], arr[7]]]
    logging.info(f"Write the info to {out}")
    with open(out, 'w') as OUT:
        print(*["Sample", "filter", "denoised", "non_overlap", "chimeric", "clean_tags"], sep="\t", file=OUT)
        for i, j in info.items():
            arr = [i, j[0] - j[1], j[1] - j[2], j[2] - j[3], j[3] - j[4], j[-1]]
            print(*arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
