#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/4/25 17:52
# @Last Modified by:   Ming
# @Last Modified time: 2023/4/25 17:52
import logging
import math
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def parse_vcf(fname):
    """
    Parse the vcf file

    :param fname: The input vcf file name
    """
    res = {}
    with open(fname, 'r') as IN:
        for line in IN:
            if not line.startswith("#"):
                arr = line.strip().split("\t")
                chromosome = arr[0]
                start = int(arr[1]) - 1
                end = arr[1]
                ref = arr[3]
                alt = arr[4]
                res[(chromosome, start, end)] = (ref, alt)
    return res


def merge_vcf(info: list, freq: float):
    """
    Merge the vcf info

    :param info: A list contain all vcf info
    :param freq: The freq of diff nucl regard as mutation
    :return res: A dict contain all mutations
    """
    num_sample = len(info) + 1
    num_mut_cutoff = math.ceil(num_sample * freq)
    logger.debug(num_mut_cutoff)
    tmp = {}
    for d in info:
        for k, v in d.items():
            if k not in tmp:
                tmp[k] = [v]
            else:
                tmp[k].append(v)

    res = {}
    for i, j in tmp.items():
        if len(j) > num_mut_cutoff:
            res[i] = j
    return res


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-d', '--din',
              required=True,
              type=click.Path(),
              help="The input dir for vcf")
@click.option('-o', '--fout',
              required=True,
              type=click.Path(),
              help="The output merge bed file")
@click.option('-f', '--freq',
              default=0.05,
              type=float,
              show_default=True,
              help="The freq of snp")
def main(din, fout, freq):
    """
    合并MicroGMT输出的结果
    """
    d_in = Path(din).absolute()
    f_out = Path(fout).absolute()

    logger.info(f"Start to read vcf file in {d_in}")
    info_vcf = []
    for f_vcf in d_in.glob("*.vcf"):
        info_vcf.append(parse_vcf(f_vcf))

    logger.info(f"Merge the vcf")
    vcf = merge_vcf(info_vcf, freq)
    logger.info(f"Output the result to {f_out}")
    with open(f_out, 'w') as OUT:
        for i, j in sorted(vcf.items(), key=lambda x: x[0][1]):
            tmp = ['|'.join(x) for x in j]
            print(*[*i, ','.join(tmp), len(tmp), "+"], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
