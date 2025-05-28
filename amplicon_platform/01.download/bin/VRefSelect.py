#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/4/19 16:58
# @Last Modified by:   Ming
# @Last Modified time: 2023/4/19 16:58
import logging
import subprocess
import sys
from pathlib import Path

import click
import numpy as np
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def get_input_stat(f_dist):
    """
    获取到每条序列距离之和最小的index

    :param f_dist: The dist file generate by mafft-distance
    """
    tmp = {}
    with open(f_dist, 'r') as IN:
        for line in IN:
            arr = line.strip().split()
            index1, index2 = arr[0].strip().split('-')
            index1 = str(int(index1) - 1)
            index2 = str(int(index2) - 1)
            dist = float(arr[1].strip().split('=')[1])
            tmp.setdefault(index1, [])
            tmp[index1].append(dist)
            tmp.setdefault(index2, [])
            tmp[index2].append(dist)
    res = {i: [np.mean(j), np.var(j)] for i, j in tmp.items()}
    return res


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--fasta',
              required=True,
              type=click.Path(),
              help="The Virus fasta genome file input")
@click.option('--out',
              default="ref.fna",
              required=True,
              type=click.Path(),
              help="The out pur ref file")
@click.option('--fstat',
              required=False,
              type=click.Path(),
              help="The out put file contain the dist info")
@click.option('--mafft',
              required=True,
              type=click.Path(),
              help="The mafft software path")
@click.option('--keep/--no-keep',
              default=False,
              show_default=True,
              help="Whether keep the middle file")
@click.option('-t', '--thread',
              default=64,
              type=int,
              show_default=True,
              help="Thread number for mafft")
def main(fasta, out, fstat, mafft, keep, thread):
    """
    从病毒基因组中挑选同其它序列距离最小的基因组作为参考序列
    """
    f_fasta = Path(fasta).absolute()
    f_out = Path(out).absolute()
    d_out = f_out.parent
    mafft = Path(mafft).absolute()

    tmp_align = d_out.joinpath("align.tmp")
    cmd_mafft = f"{mafft} --auto --quiet --thread {thread} {f_fasta} > {tmp_align}"
    logger.info(f"Start to do mafft align")
    try:
        (mafft_out, mafft_error) = subprocess.Popen(cmd_mafft, shell=True, stdout=subprocess.PIPE).communicate()
    except:
        sys.exit(logger.error(mafft_error))

    tmp_dist = d_out.joinpath("dist.tmp")
    cmd_dist = f"{mafft}-distance {tmp_align} > {tmp_dist}"
    logger.info(f"Start to get the distance")
    try:
        (dist_out, dist_error) = subprocess.Popen(cmd_dist, shell=True, stdout=subprocess.PIPE).communicate()
    except:
        sys.exit(logger.error(dist_error))

    logger.info(f"Get the reference to file: {f_out}")
    info_dist = get_input_stat(tmp_dist)
    ref_index = int(sorted(info_dist.items(), key=lambda x: x[1])[0][0])
    res_dist = {}
    with open(f_out, 'w') as OUT:
        index = 0
        for record in SeqIO.parse(f_fasta, 'fasta'):
            res_dist[record.description] = info_dist[str(index)]
            if index == ref_index:
                print(f">{record.description}\n{record.seq}", file=OUT)
            index += 1

    if fstat:
        logger.info(f"Out put the stat info to file: {fstat}")
        with open(fstat, 'w') as OUT:
            for i, j in res_dist.items():
                print(*[i, *j], sep="\t", file=OUT)

    if keep == False:
        tmp_align.unlink(missing_ok=True)
        tmp_dist.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
