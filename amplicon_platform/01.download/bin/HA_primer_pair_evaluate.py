#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/26 15:30
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/26 15:30
import logging
import re
from itertools import product
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def is_pair(f_info, r_info):
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
        logger.warning("Multi amplify find")
        return False


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--primer_name',
              required=True,
              type=click.Path(),
              help="The primer name info file")
@click.option('--blast',
              required=True,
              type=click.Path(),
              help="The custom blast result")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The out put result")
@click.option('--info',
              required=True,
              type=click.Path(),
              help="The info file contain the HA segment number in database")
def main(primer_name, blast, out, info):
    """
    计算成对引物的特异性
    """
    f_primer_name = Path(primer_name).absolute()
    f_blast = Path(blast).absolute()
    f_out = Path(out).absolute()
    f_info = Path(info).absolute()

    logger.info(f"Parse the blast result")
    # TODO: 这里比较消耗内存，后续优化
    info_blast = {}
    with open(f_blast, 'r') as IN:
        for line in IN:
            arr = line.strip().split('\t')
            info_blast.setdefault(arr[0], {})
            info_blast[arr[0]].setdefault(arr[1], [])
            info_blast[arr[0]][arr[1]].append(arr[2:])

    logger.info("Parse the HA Segment Number info")
    target = ["H1", "H3", "H5", "H7", "H9"]
    info_num = {i: 0 for i in target}
    with open(f_info, 'r') as IN:
        for line in IN:
            arr = line.strip().split("\t")
            if arr[0] in info_num:
                info_num[arr[0]] = int(arr[1])

    logger.info(f"Start to analysis and output the result to {f_out}")
    res = {}
    with open(f_primer_name, 'r') as IN, open(f_out, 'w') as OUT:
        print(*["Name", "F", "R", "H1 Number(percent %)", "H3 Number(percent %)", "H5 Number(percent %)",
                "H7 Number(percent %)", "H9 Number(percent %)"], sep="\t", file=OUT)
        for line in IN:
            arr = line.strip().split("\t")
            name = arr[0]
            res[name] = {i: 0 for i in target}
            primer_f, primer_r = arr[1], arr[2]
            for subject, subject_f_info in info_blast[primer_f].items():
                ha_type = re.findall(r"A_/_(H\d+)", subject)[0]
                if subject in info_blast[primer_r]:
                    subject_r_info = info_blast[primer_r][subject]
                    if len(subject_f_info) > 1:
                        logger.warning(f"{name} primer F: {primer_f} has align multi place on {subject}")
                        logger.debug(subject_f_info)
                    if len(subject_r_info) > 1:
                        logger.warning(f"{name} primer R: {primer_r} has align multi place on {subject}")
                        logger.debug(subject_r_info)
                    if is_pair(subject_f_info, subject_r_info):
                        if ha_type in info_num:
                            res[name][ha_type] += 1
            content = [name, primer_f, primer_r]
            for ha_type in target:
                num = res[name][ha_type]
                percent = round(num / info_num[ha_type] * 100, 2) if info_num[ha_type] != 0 else "0.00"
                content.append(f"{num}({percent})")
            print(*content, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
