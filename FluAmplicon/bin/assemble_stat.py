#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/24 9:35
# @Last Modified by:   Ming
# @Last Modified time: 2022/6/24 9:35
import os

from Bio import SeqIO
import logging
import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-d', '--d_in',
              required=True,
              type=click.Path(),
              help="The input dir contain segment fasta file")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output stat file")
def main(d_in, out):
    """
    统计IRMA组装结果内8个片段的信息
    """
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    res = []
    for segment in segments:
        f_in = os.path.join(d_in, f"{segment}.fasta")
        if os.path.exists(f_in):
            logging.info(f"Read {segment} file {f_in}")
            for record in SeqIO.parse(f_in, "fasta"):
                res.append(len(record.seq))
        else:
            res.append(0)
    with open(out, 'w') as OUT:
        print(*segments, sep="\t", file=OUT)
        print(*res, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
