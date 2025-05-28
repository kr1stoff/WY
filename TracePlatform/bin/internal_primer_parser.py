#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/2/10 16:46
# @Last Modified by:   mengxf
# @Last Modified time: 2023/5/11
import logging
import sys
from pathlib import Path

import click
from Bio import SeqIO

d_bin = Path(__file__).absolute().parent
d_lib = d_bin.parent.joinpath("lib")

sys.path.append(str(d_lib))
from Primer3 import Parser

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-f', '--file',
              required=True,
              type=click.Path(),
              help="The Primer3 out put internal primer result")
@click.option('--ref',
              required=True,
              type=click.Path(),
              help="The reference file")
@click.option('--left',
              default=50,
              show_default=True,
              type=int,
              help="The left distance to the start")
@click.option('--right',
              default=50,
              show_default=True,
              type=int,
              help="The right distance to the end")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
@click.option('--opt_prm_len',
              default=20,
              type=int,
              show_default=True,
              help="The optimal primer length.")
@click.option('--len_mult',
              default=1.0,
              type=float,
              show_default=True,
              help="The primer length penalty multiplier.")
@click.option('--min_amp_len',
              default=70,
              type=int,
              show_default=True,
              help="The minimum amplicon length.")
@click.option('--min_amp_mult',
              default=5,
              type=int,
              show_default=True,
              help="The primer length penalty multiplier.")
def main(file, ref, left, right, out,
         opt_prm_len, len_mult, min_amp_len, min_amp_mult):
    """
    解析使用Primer3设计的内引物结果文件, 使用PrimerExpress罚分规则.
    """
    f_primer3 = Path(file).absolute()
    f_ref = Path(ref).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True)
    logging.info(f"Parse the reference file")
    for record in SeqIO.parse(f_ref, "fasta"):
        ref_length = len(record.seq)
    logging.info(f"Parse the file {f_primer3}")
    res = Parser("internal", opt_prm_len, len_mult, min_amp_len, min_amp_mult)
    res.load_file(f_primer3)
    # res.filter()
    res.filter_region(left, ref_length - right)
    # 引物信息列表文件
    f_tsv = d_out.joinpath("primer.xls")
    # 引物及探针的序列
    f_fasta = d_out.joinpath("primer.fasta")
    # bed 文件
    f_bed = d_out.joinpath("primer.bed")
    res.to_csv(f_tsv, sep="\t")
    res.to_fasta(f_fasta)
    res.to_bed(f_bed)


if __name__ == "__main__":
    main()
