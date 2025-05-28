#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/8/7 19:46
# @Last Modified by:   Ming
# @Last Modified time: 2023/8/7 19:46
import hashlib
import logging
from itertools import product
from pathlib import Path

import click
from Bio import Seq

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def get_md5(sequence: str):
    """
    获取字符串的md5值
    """
    res = hashlib.md5(sequence.encode("utf-8")).hexdigest()
    return res


def extend_ambiguous_dna(seq):
    """return list of all possible sequences given an ambiguous DNA input"""
    d = Seq.IUPACData.ambiguous_dna_values
    return list(map("".join, product(*map(d.get, seq))))


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--data',
              required=True,
              type=click.Path(),
              help="The input sequence primer.xls file or the file after get the degenerate")
@click.option('-o', "--out",
              required=True,
              type=click.Path(),
              help="The output dir")
def main(data, out):
    """
    使用引物设计的结果文件生成后续评估所需要的一系列文件
    """
    f_data = Path(data).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)

    f_sequence = d_out.joinpath("sequence.fasta")
    f_primer_info = d_out.joinpath("primer.info.tsv")

    logger.info(f"Parse the input data file {f_data} and output the info")
    pool = set()
    with open(f_data, 'r') as IN, open(f_sequence, 'w') as SEQ, open(f_primer_info, 'w') as INTERNAL:
        next(IN)
        for line in IN:
            arr = line.strip().split("\t")
            index = arr[0]
            primer_f = arr[7].upper()
            primer_r = arr[13].upper()
            # 引物/探针名称
            primer_f_md5 = get_md5(primer_f)
            primer_r_md5 = get_md5(primer_r)

            seqs = [primer_f, primer_r]
            md5s = [primer_f_md5, primer_r_md5]
            for i, j in zip(md5s, seqs):
                if i not in pool:
                    for split_seq in extend_ambiguous_dna(j):
                        print(f">{i}\n{split_seq}", file=SEQ)
                    pool.add(i)
            print(*[index, primer_f_md5, primer_r_md5], sep="\t", file=INTERNAL)


if __name__ == "__main__":
    main()
