#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/2 9:57
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/2 9:57
import logging
from pathlib import Path

import click
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def get_gene_name(title):
    """
    Get the gene name of the segment
    """
    res = title.strip().split("|")[2]
    return res


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-f', '--fasta',
              required=True,
              type=click.Path(),
              help="The input fasta file from GISAID flu database")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put dir")
def main(fasta, out):
    """
    Split fasta file download from GISAID Flu database to split genes

    The gene name should use the format below
    >1165324|B/Minnesota/01/2018|NP|Original|B_/_H0N0
    attttcttgtgaacttcaagcaccagtaaaagaactgaaaatcaaaatgtccaacatgg
    """
    f_fasta = Path(fasta).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)

    genes = ["NP", "NS", "MP", "PA", "PB2", "PB1", "NA", "HA"]
    info = {i: {} for i in genes}
    logger.info(f"Parse the input fasta file: {f_fasta}")
    for record in SeqIO.parse(f_fasta, "fasta"):
        name = record.name
        seq = record.seq
        _type = get_gene_name(name)
        info[_type][name] = seq

    logger.info(f"Output the result info to dir: {d_out}")
    for i, j in info.items():
        f_out = d_out.joinpath(f"{i}.fasta")
        with open(f_out, 'w') as OUT:
            for name, seq in j.items():
                print(f">{name}\n{seq}", file=OUT)


if __name__ == "__main__":
    main()
