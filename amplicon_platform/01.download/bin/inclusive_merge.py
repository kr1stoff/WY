#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/19 15:43
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/19 15:43
import hashlib
import logging
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def parse_single(f_name):
    """
    Parse the probe inclusive result file
    """
    res = {}
    with open(f_name, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split("\t")
            res[arr[0]] = arr[2]

    return res


def parse_pair(f_name):
    """
    Parse the primer inclusive result file
    """
    res = {}
    with open(f_name, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split("\t")
            res[arr[0]] = arr[3:]
    return res


def get_md5(sequence: str):
    """
    获取字符串的md5值
    """
    res = hashlib.md5(sequence.encode("utf-8")).hexdigest()
    return res


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.1")
@click.option('--probe',
              required=True,
              type=click.Path(),
              help="The probe specificity result file")
@click.option('--internal',
              required=True,
              type=click.Path(),
              help="The internal specificity result file")
@click.option('--external',
              required=True,
              type=click.Path(),
              help="The external specificity result file")
@click.option('--info',
              required=True,
              type=click.Path(),
              help="The pPCR pipe result or after get the degenerate")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put file name")
def main(probe, internal, external, info, out):
    """
    合并包容性分析的结果

    **整数部分为可以成功扩增的数目，小数为占数据库基因组的百分比**
    """
    f_probe = Path(probe).absolute()
    f_internal = Path(internal).absolute()
    f_external = Path(external).absolute()
    f_info = Path(info).absolute()
    f_out = Path(out).absolute()

    logger.info(f"Parse the probe specificity result: {f_probe}")
    info_probe = parse_single(f_probe)
    logger.debug(info_probe)
    logger.info(f"Parse the internal primer specificity result: {f_internal}")
    info_internal = parse_pair(f_internal)
    logger.debug(info_internal)
    logger.info(f"Parse the external specificity result: {f_external}")
    info_external = parse_pair(f_external)
    logger.debug(info_external)

    logger.info(f"Parse the info file {f_info} and output the result to {f_out}")
    with open(f_info, 'r') as IN, open(f_out, 'w') as OUT:
        next(IN)
        header = ["Index", "InternaleF(%)", "InternalR(%)", "InternalPair(%)", "Probe(%)", "ExternalF(%)",
                  "ExternalR(%)", "ExternalPair(%)"]
        print(*header, sep="\t", file=OUT)
        for line in IN:
            arr = line.strip().split("\t")
            index = arr[0]
            internal_content = info_internal[index] if index in info_internal else ["0(0.00)", "0(0.00)", "0(0.00)"]
            probe_content = info_probe[index] if index in info_probe else "0(0.00)"
            external_content = info_external[index] if index in info_external else ["0(0.00)", "0(0.00)", "0(0.00)"]
            new_arr = [index, *internal_content, probe_content, *external_content]
            print(*new_arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
