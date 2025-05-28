#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/14 14:32
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/14 14:32
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
@click.option('--dirin',
              required=True,
              type=click.Path(),
              help="The input dir")
@click.option('-p', '--pattern',
              required=True,
              default="_genomic.fna",
              show_default=True,
              type=click.Path(),
              help="The pattern to use for input genome")
@click.option('-t', '--taxon',
              required=True,
              type=click.Path(),
              help="The taxon id")
@click.option('-o', '--out',
              required=True,
              default="assembly_chromosome.tsv",
              show_default=True,
              type=click.Path(),
              help="The out put file name")
def main(dirin, pattern, taxon, out):
    """
    生成数据库所需的assembly_chromosome.tsv文件(只针对从refseq和Genbank下载的数据)
    """
    d_in = Path(dirin).absolute()
    f_out = Path(out).absolute()

    logger.info(f"You input dir is: {d_in}")
    logger.info(f"You pattern is: {pattern}")
    logger.info(f"You taxonid is: {taxon}")

    header = ["taxon_id", "genome_name", "chromosome_included"]
    with open(f_out, 'w') as OUT:
        print(*header, sep="\t", file=OUT)
        for i in d_in.glob(f"*{pattern}"):
            name = i.name
            chromosomes = []
            for record in SeqIO.parse(i, "fasta"):
                chromosomes.append(record.name)
            print(*[taxon, name, ','.join(chromosomes)], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
