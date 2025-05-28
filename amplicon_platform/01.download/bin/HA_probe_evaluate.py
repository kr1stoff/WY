#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/26 17:47
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/26 17:47
import logging
import re
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--probe_name',
              required=True,
              type=click.Path(),
              help="The probe name info file")
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
def main(probe_name, blast, out, info):
    """
    计算探针的特异性
    """
    f_probe_name = Path(probe_name).absolute()
    f_blast = Path(blast).absolute()
    f_out = Path(out).absolute()
    f_info = Path(info).absolute()

    logger.info(f"Parse the probe name")
    name2probe = {}
    pool_probe = set()
    with open(f_probe_name, 'r') as IN:
        for line in IN:
            arr = line.strip().split('\t')
            name2probe[arr[0]] = arr[1]
            pool_probe.add(arr[1])

    logger.info("Parse the HA Segment Number info")
    target = ["H1", "H3", "H5", "H7", "H9"]
    info_num = {i: 0 for i in target}
    with open(f_info, 'r') as IN:
        for line in IN:
            arr = line.strip().split("\t")
            if arr[0] in info_num:
                info_num[arr[0]] = int(arr[1])

    logger.info(f"Parse the blast result")
    info_blast = {}
    with open(f_blast, 'r') as IN:
        for line in IN:
            arr = line.strip().split('\t')
            query = arr[0]
            if query in pool_probe:
                subject = arr[1]
                ha_type = re.findall(r"A_/_(H\d+)", subject)[0]
                info_blast.setdefault(query, {i: 0 for i in target})
                if ha_type in info_blast[query]:
                    info_blast[query][ha_type] += 1

    logger.debug(info_blast)

    logger.info(f"Output the result to {f_out}")
    with open(f_out, 'w') as OUT:
        print(*["Name", "Probe", "H1 Number(percent %)", "H3 Number(percent %)", "H5 Number(percent %)",
                "H7 Number(percent %)", "H9 Number(percent %)"], sep="\t", file=OUT)
        for i, j in name2probe.items():
            logger.debug(f"({i},{j})")
            content = [i, j]
            for ha_type in target:
                num = info_blast[j][ha_type] if j in info_blast and ha_type in info_blast[j] else 0
                percent = round(num / info_num[ha_type] * 100, 2) if info_num[ha_type] != 0 else "0.00"
                content.append(f"{num}({percent})")
            print(*content, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
