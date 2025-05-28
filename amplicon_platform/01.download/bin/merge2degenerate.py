#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/14 14:38
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/14 14:38
import logging
import sys
from pathlib import Path

import click
from Bio import Seq

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
#### Some Global Variable
# 兼并碱基列表
dict_degenerate = {('A',): 'A',
                   ('C',): 'C',
                   ('G',): 'G',
                   ('T',): 'T',
                   ('A', 'G'): 'R',
                   ('C', 'T'): 'Y',
                   ('C', 'G'): 'S',
                   ('A', 'T'): 'W',
                   ('G', 'T'): 'K',
                   ('A', 'C'): 'M',
                   ('C', 'G', 'T'): 'B',
                   ('A', 'G', 'T'): 'D',
                   ('A', 'C', 'T'): 'H',
                   ('A', 'C', 'G'): 'V',
                   ('A', 'C', 'G', 'T'): 'N'}


#### Some Function
def get_degenerate(l_alpha: list):
    """
    获取对应碱基列表对应的兼并碱基

    :param l_alpha: A list contain the A、T、G、C
    """
    key = tuple(sorted([i.upper() for i in l_alpha]))
    if key in dict_degenerate:
        return dict_degenerate[tuple(sorted(l_alpha))]
    else:
        sys.exit(logger.error(f"{key} not in the degenerate table, please check it"))


def seq2degenerate(seq: str, chromosome, start, info: dict, strand: str = "+"):
    """
    将提供的序列内的碱基转变未兼并碱基

    :param seq: The primer sequence
    :param chromosome: The chromosome name
    :param start: The start position of the primer(0 base)
    :param info: The info dict contain the snp info
    :param strand: The primer is in the + strand or the - strand
    """
    start = int(start)
    base = []
    if strand == '+':
        pass
    elif strand == '-':
        seq = Seq.Seq(seq).reverse_complement()
    else:
        sys.exit(logger.error(f"strand param should only be '+' or '-', but your input is {strand}"))

    for i in range(len(seq)):
        key = tuple([chromosome, str(start + i)])
        if key in info:
            base.append(get_degenerate(info[key]))
        else:
            base.append(seq[i])

    res = ''.join(base)
    if strand == '-':
        res = str(Seq.Seq(res).reverse_complement())

    return res


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--merge',
              required=True,
              type=click.Path(),
              help="The merge.xls file generate by pCPR pipe")
@click.option('--snp',
              required=True,
              type=click.Path(),
              help="The all.variants.bed generate by the pPCR pipe")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output file")
@click.option('--freq',
              default=0.05,
              show_default=True,
              type=float,
              help="The min freq of a site consider as snp")
def main(merge, snp, out, freq):
    """
    将引物设计结果的merge.xls文件内的碱基按照条件替换为兼并碱基

    *注意流程中的突变文件是有一定的筛选条件的，如果freq低于流程中的筛选条件，可按照freq=0的条件重新生成all.variants.bed文件*
    """
    f_merge = Path(merge).absolute()
    f_snp = Path(snp).absolute()
    f_out = Path(out).absolute()

    logger.info(f"Parse the SNP info from file {f_snp}")
    info_snp = {}
    with open(f_snp, 'r') as IN:
        for line in IN:
            arr = line.strip().split("\t")
            key = tuple([arr[0], arr[1]])
            num_total = int(arr[5])
            for type_, num_ in zip(arr[3].split(','), arr[4].split(',')):
                num_ = int(num_)
                if num_ / num_total > freq:
                    info_snp.setdefault(key, set())
                    for i in type_.strip().split('|'):
                        info_snp[key].add(i)

    logger.info(f"Parse the pPCR result file {f_merge}")
    logger.info(f"The result will out put to file {f_out}")
    with open(f_merge, 'r') as IN, open(f_out, 'w') as OUT:
        print(next(IN).strip(), file=OUT)
        for line in IN:
            arr = line.strip().split('\t')
            chromosome = arr[2]
            # 内引物F
            start_forward_internal = arr[3]
            seq_forward_internal = arr[8]
            arr[8] = seq2degenerate(seq_forward_internal, chromosome, start_forward_internal, info_snp)
            # 内引物R
            start_reverse_internal = arr[9]
            seq_reverse_internal = arr[14]
            arr[14] = seq2degenerate(seq_reverse_internal, chromosome, start_reverse_internal, info_snp, strand='-')
            # 探针
            start_probe = arr[15]
            seq_probe = arr[20]
            arr[20] = seq2degenerate(seq_probe, chromosome, start_probe, info_snp)
            # 外引物F
            start_forward_external = arr[24]
            seq_forward_external = arr[29]
            arr[29] = seq2degenerate(seq_forward_external, chromosome, start_forward_external, info_snp)
            # 外引物R
            start_reverse_external = arr[30]
            seq_reverse_external = arr[35]
            arr[35] = seq2degenerate(seq_reverse_external, chromosome, start_reverse_external, info_snp, strand='-')
            print(*arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
