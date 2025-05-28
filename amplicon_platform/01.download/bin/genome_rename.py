#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/27 19:19
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/27 19:19
import logging
from pathlib import Path

import click
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--fasta',
              required=True,
              type=click.Path(),
              help="The input fasta file")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The output fasta file")
def main(fasta, out):
    """
    将fasta内的序列名称从0开始命名，从而确保SNP时不会出现错误
    """
    f_fasta = Path(fasta).absolute()
    f_out = Path(out).absolute()

    with open(f_out, 'w') as OUT:
        flag = 0
        for record in SeqIO.parse(f_fasta, "fasta"):
            print(f">{flag}\n{record.seq}", file=OUT)
            flag += 1


if __name__ == "__main__":
    main()
