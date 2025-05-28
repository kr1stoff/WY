#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/14 14:58
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/14 14:58
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
@click.option('-f', '--fasta',
              required=True,
              type=click.Path(),
              help="The input fasta genome file from nucleotide")
@click.option('-l', '--glength',
              required=True,
              type=int,
              help="The represent genome length")
@click.option('-r', '--grange',
              default=10,
              show_default=True,
              type=int,
              help="The genome range with represent genome length(+-)")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put file name")
def main(fasta, glength, grange, out):
    """
    对从nucleotide下载的病毒序列按照指定长度的上下浮动范围进行过滤
    """
    f_in = Path(fasta).absolute()
    f_out = Path(out).absolute()

    l_min = glength * (1 - grange / 100)
    l_max = glength * (1 + grange / 100)
    logger.info(f"Your target genome length is between {l_min} ~ {l_max}")
    with open(f_out, 'w') as OUT:
        for record in SeqIO.parse(f_in, "fasta"):
            if l_min < len(record.seq) < l_max:
                print(f">{record.description}\n{record.seq}", file=OUT)


if __name__ == "__main__":
    main()
