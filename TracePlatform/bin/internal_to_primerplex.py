#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/9 9:24
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/13 9:24
import logging

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-i', '--finput',
              required=True,
              type=click.Path(),
              help="The bed file of primer")
@click.option('-o', '--foutput',
              required=True,
              type=click.Path(),
              help="The bed file of primerplex")
def main(finput, foutput):
    """
    将引物的bed文件生产primerplex的格式
    """
    logger.info(f"Get info from {finput}")
    with open(finput, 'r') as IN, open(foutput, 'w') as OUT:
        for line in IN:
            arr1 = line.strip().split("\t")
            arr2 = next(IN).strip().split('\t')
            next(IN)
            res = [arr1[0], arr1[1], arr2[2], arr1[3].strip().rsplit('-', 1)[0],
                   1, 'B', 'W']
            print(*res, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
