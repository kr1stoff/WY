#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/7/6 11:57
# @Last Modified by:   Ming
# @Last Modified time: 2022/7/6 11:57
from Bio import SeqIO
from pathlib import Path
import logging
import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-s', '--scheme',
              required=True,
              type=click.Path(),
              help="The mlst scheme file in the species dir")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put fasta file")
def main(scheme, out):
    """
    使用 mlst 软件的 pubmlst 文件夹内某个物种的序列生成 blast 目录所需的序列
    """
    scheme = Path(scheme)
    d_scheme = scheme.parent
    name_scheme = scheme.stem
    logging.info(f"Parse the file: {scheme}")
    genes = []
    with open(scheme, 'r') as IN:
        genes = next(IN).strip().split("\t")[1:]

    with open(out, 'w') as OUT:
        for gene in genes:
            f_gene = d_scheme.joinpath(f"{gene}.tfa")
            if f_gene.exists():
                for record in SeqIO.parse(f_gene, "fasta"):
                    name = f"{name_scheme}.{record.name}"
                    print(f">{name}\n{record.seq}", file=OUT)


if __name__ == "__main__":
    main()
