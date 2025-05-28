#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/17 16:59
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/17 16:59
import logging
from itertools import product
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def is_pair(f_info, r_info, multi=False):
    """
    判断F引物同R引物能否配对
    """
    num = 0
    logger.debug(f_info)
    logger.debug(r_info)
    for f, r in product(f_info, r_info):
        f_start, f_end, f_strand = int(f[6]), int(f[7]), f[-1]
        r_start, r_end, r_strand = int(r[6]), int(r[7]), r[-1]
        if f_strand != r_strand:
            if f_strand == "plus":
                len_amp = r_start - f_start + 1
            else:
                len_amp = f_start - r_start + 1
            if 0 < len_amp < 2000:
                num += 1
    logger.debug(num)
    if num == 0:
        return False
    elif num == 1:
        return True
    else:
        if multi:
            return True
        else:
            logger.warning("Multi amplify find")
            return False


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--blast',
              required=True,
              type=click.Path(),
              help="The blast result for primer/probe")
@click.option('--info',
              required=True,
              type=click.Path(),
              help="The primer/probe info file")
@click.option('-t', '--type_',
              default="primer",
              show_default=True,
              type=click.Choice(["primer", "probe"]),
              help="The input type")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output file name")
@click.option('--chromosome_info',
              required=True,
              type=click.Path(),
              help="The assembly_chromosome.tsv file to use")
@click.option('--multi',
              is_flag=True,
              default=False,
              show_default=True,
              help="Whether allow multi alignment on one chromosome(useful for degenerate)")
def main(blast, info, type_, out, chromosome_info, multi):
    """
    引物包容性评估

    **此脚本并不会对blast结果进行判定，故须提前对blast结果进行过滤**
    **如果有blast结果为兼并碱基的展开结果，请加入--multi参数**
    """
    f_blast = Path(blast).absolute()
    f_info = Path(info).absolute()
    f_chromosome = Path(chromosome_info).absolute()
    f_out = Path(out).absolute()

    logger.info(f"Parse the blast result")
    # TODO: 这里比较消耗内存，后续优化
    info_blast = {}
    with open(f_blast, 'r') as IN:
        for line in IN:
            arr = line.strip().split('\t')
            info_blast.setdefault(arr[0], {})
            info_blast[arr[0]].setdefault(arr[1], [])
            info_blast[arr[0]][arr[1]].append(arr[2:])

    total_genome_num = 0
    chromosome2genome = {}
    logger.info(f"Parse the chromosome info file: {f_chromosome}")
    with open(f_chromosome, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split("\t")
            for i in arr[2].strip().split(','):
                chromosome2genome[i] = arr[1]
            total_genome_num += 1

    logger.info(f"Start to analysis and output the result to {f_out}")
    with open(f_info, 'r') as IN, open(f_out, 'w') as OUT:
        if type_ == "probe":
            header = ["Name", "Probe", "Inclusive Number(percent %)"]
            print(*header, sep="\t", file=OUT)
            for line in IN:
                arr = line.strip().split("\t")
                name = arr[0]
                probe = arr[1]
                pool = set()
                if probe in info_blast:
                    for i in info_blast[probe].keys():
                        pool.add(chromosome2genome[i])
                percent = len(pool) / total_genome_num * 100
                print(*[name, probe, f"{len(pool)}({percent:.2f})"], sep="\t", file=OUT)
        else:
            header = ["Name", "Primer F", "Primer R", "Inclusive F Number(percent %)", "Inclusive R Number(percent %)",
                      "Inclusive Pair Number(percent %)"]
            print(*header, sep="\t", file=OUT)
            for line in IN:
                arr = line.strip().split("\t")
                name = arr[0]
                primer_f, primer_r = arr[1], arr[2]
                pool_f = set(
                    [chromosome2genome[i] for i in info_blast[primer_f].keys()]) if primer_f in info_blast else set()
                pool_r = set(
                    [chromosome2genome[i] for i in info_blast[primer_r].keys()]) if primer_r in info_blast else set()
                pool_pair = set()
                if primer_f in info_blast and primer_r in info_blast:
                    for subject, subject_f_info in info_blast[primer_f].items():
                        if subject in info_blast[primer_r]:
                            subject_r_info = info_blast[primer_r][subject]
                            if not multi:
                                if len(subject_f_info) > 1:
                                    logger.warning(f"{name} primer F: {primer_f} has align multi place on {subject}")
                                    logger.debug(subject_f_info)
                                if len(subject_r_info) > 1:
                                    logger.warning(f"{name} primer R: {primer_r} has align multi place on {subject}")
                                    logger.debug(subject_r_info)
                            if is_pair(subject_f_info, subject_r_info, multi):
                                pool_pair.add(chromosome2genome[subject])
                num_f = len(pool_f)
                percent_f = num_f / total_genome_num * 100
                num_r = len(pool_r)
                percent_r = num_r / total_genome_num * 100
                num_pair = len(pool_pair)
                percent_pair = num_pair / total_genome_num * 100
                print(*[name, primer_f, primer_r, f"{num_f}({percent_f:.2f})", f"{num_r}({percent_r:.2f})",
                        f"{num_pair}({percent_pair:.2f})"], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
