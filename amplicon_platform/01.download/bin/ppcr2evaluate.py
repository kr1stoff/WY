#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/28 10:04
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/28 10:04
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
@click.version_option(version="v1.1.0")
@click.option('--data',
              required=True,
              type=click.Path(),
              help="The input pPCR primer.xls file or the file after get the degenerate")
@click.option('-o', "--out",
              required=True,
              type=click.Path(),
              help="The output dir")
def main(data, out):
    """
    使用pPCR的结果文件生成后续评估所需要的一系列文件
    """
    f_data = Path(data).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)

    f_sequence = d_out.joinpath("sequence.fasta")
    f_internal_info = d_out.joinpath("internal.info.tsv")
    f_probe_info = d_out.joinpath("probe.info.tsv")
    f_external_info = d_out.joinpath("external.info.tsv")

    logger.info(f"Parse the input data file {f_data} and output the info")
    pool = set()
    with open(f_data, 'r') as IN, open(f_sequence, 'w') as SEQ, open(f_internal_info, 'w') as INTERNAL, open(
            f_probe_info, 'w') as PROBE, open(f_external_info, 'w') as EXTERNAL:
        next(IN)
        for line in IN:
            arr = line.strip().split("\t")
            index = arr[0]
            primer_internal_f = arr[8].upper()
            primer_internal_r = arr[14].upper()
            probe = arr[20].upper()
            primer_external_f = arr[29].upper()
            primer_external_r = arr[35].upper()
            # 引物/探针名称
            primer_internal_f_md5 = get_md5(primer_internal_f)
            primer_internal_r_md5 = get_md5(primer_internal_r)
            probe_md5 = get_md5(probe)
            primer_external_f_md5 = get_md5(primer_external_f)
            primer_external_r_md5 = get_md5(primer_external_r)

            seqs = [primer_internal_f, primer_internal_r, probe, primer_external_f, primer_external_r]
            md5s = [primer_internal_f_md5, primer_internal_r_md5, probe_md5, primer_external_f_md5,
                    primer_external_r_md5]
            for i, j in zip(md5s, seqs):
                if i not in pool:
                    for split_seq in extend_ambiguous_dna(j):
                        print(f">{i}\n{split_seq}", file=SEQ)
                    pool.add(i)
            print(*[index, primer_internal_f_md5, primer_internal_r_md5], sep="\t", file=INTERNAL)
            print(*[index, probe_md5], sep="\t", file=PROBE)
            print(*[index, primer_external_f_md5, primer_external_r_md5], sep="\t", file=EXTERNAL)


if __name__ == "__main__":
    main()
