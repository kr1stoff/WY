#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/15 16:42
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/15 16:42
import logging
from pathlib import Path

import click
import numpy as np
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-m', '--msa',
              required=True,
              type=click.Path(),
              help="The fasta format MSA result")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put result")
def main(msa, out):
    """
    对多序列比对的结果进行打分，从而指导参考基因组的选择
    """
    f_msa = Path(msa).absolute()
    f_out = Path(out).absolute()

    names = []
    seqs = []
    logger.info(f"Read the MSA file: {f_msa}")
    for record in SeqIO.parse(f_msa, "fasta"):
        name = record.name
        seq = record.seq
        names.append(name)
        seqs.append(np.array(list(seq)))
    seqs = np.array(seqs)

    logger.info("Get the score")
    nrow, ncol = seqs.shape
    res = np.zeros((nrow, ncol))
    for i in range(ncol):
        col = seqs[:, i]
        unique, counts = np.unique(col, return_counts=True)
        ratio = counts / len(col)
        dic = dict(zip(unique, ratio))
        for j in range(nrow):
            res[j, i] = dic[seqs[j, i]]

    logger.info(f"Output the result to: {f_out}")
    with open(f_out, 'w') as OUT:
        for i, j in sorted(zip(names, res.sum(axis=1)), key=lambda x: x[1], reverse=True):
            print(*[i, j], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
