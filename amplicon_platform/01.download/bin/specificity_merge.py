#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/19 11:19
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/19 11:19
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
    Parse the probe specificity result file
    """
    res = {}
    with open(f_name, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split("\t")
            res[arr[0]] = f"{arr[2]}({arr[1]})"

    return res


def parse_pair(f_name):
    """
    Parse the primer specificity result file
    """
    res = {}
    with open(f_name, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split("\t")
            res[arr[0]] = [f"{arr[4]}({arr[3]})", f"{arr[6]}({arr[5]})", f"{arr[2]}({arr[1]})"]
    return res


def get_md5(sequence: str):
    """
    获取字符串的md5值
    """
    res = hashlib.md5(sequence.encode("utf-8")).hexdigest()
    return res


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
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
    合并特异性分析的结果文件

    整合后的结果中整数部分表示获得的比对的结果数，百分比表示其中有多少条比对到目标物种
    """
    f_probe = Path(probe).absolute()
    f_internal = Path(internal).absolute()
    f_external = Path(external).absolute()
    f_info = Path(info).absolute()
    f_out = Path(out).absolute()

    logger.info(f"Parse the probe specificity result: {f_probe}")
    info_probe = parse_single(f_probe)
    logger.info(f"Parse the internal primer specificity result: {f_internal}")
    info_internal = parse_pair(f_internal)
    logger.debug(info_internal)
    logger.info(f"Parse the external specificity result: {f_external}")
    info_external = parse_pair(f_external)

    logger.info(f"Parse the info file {f_info} and output the result to {f_out}")
    with open(f_info, 'r') as IN, open(f_out, 'w') as OUT:
        next(IN)
        header = ["Index", "InternaleF(%)", "InternalR(%)", "InternalPair(%)", "Probe(%)", "ExternalF(%)",
                  "ExternalR(%)", "ExternalPair(%)"]
        print(*header, sep="\t", file=OUT)
        for line in IN:
            arr = line.strip().split("\t")
            index = arr[0]
            logger.debug(arr[20])
            probe_md5 = get_md5(arr[20].upper())
            probe_content = info_probe[probe_md5] if probe_md5 in info_probe else "0(0.00)"
            new_arr = [index, *info_internal[index], probe_content, *info_external[index]]
            print(*new_arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
