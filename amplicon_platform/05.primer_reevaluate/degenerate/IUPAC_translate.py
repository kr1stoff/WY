#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/10 15:22
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/10 15:22
import logging
from itertools import product

import click
from Bio import Seq, SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


# Some Function
def extend_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = Seq.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


def primers_to_fasta(name, seq_list, add):
    """return fasta string of primers with tracing newline"""
    fas = ""
    for i in range(len(seq_list)):
        if (add):
            fas += f">{name}[{i}]\n{seq_list[i]}\n"  # 展开后id后面加[1],[2]等
        else:
            fas += f">{name}\n{seq_list[i]}\n"  # 展开后id后面不加[1],[2]等
    return fas


# Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fastain',
              required=True,
              type=click.Path(),
              help="The input primer fasta file")
@click.option('-o', '--fastaout',
              required=True,
              type=click.Path(),
              help="The output primer fasta file")
@click.option('--add/--no-add',
              default=False,
              show_default=True,
              help="是否在ID后面加副ID")
def cli(fastain, fastaout, add):
    """
    将带有兼并碱基的序列展开为ATGC(展开后的序列带有相同的名称)
    """
    with open(fastain, "r") as fin, open(fastaout, "w") as fout:
        for record in SeqIO.parse(fin, "fasta"):
            explicit = extend_ambiguous_dna(record.seq.upper())
            fasta = primers_to_fasta(record.id, explicit, add)
            fout.write(fasta)


if __name__ == "__main__":
    cli()
