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
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put dir")
def main(fasta, out):
    """
    对从nucleotide下载的病毒分割到单独的文件
    """
    f_in = Path(fasta).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True)

    logger.info(f"Split your sequence to dir: {d_out}")
    for record in SeqIO.parse(f_in, "fasta"):
        name = record.name
        f_out = d_out.joinpath(f"{name}_genomic.fna")
        with open(f_out, 'w') as OUT:
            print(f">{record.description}\n{record.seq}", file=OUT)


if __name__ == "__main__":
    main()
