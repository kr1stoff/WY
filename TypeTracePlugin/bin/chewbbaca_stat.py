#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/7/11 13:55
# @Last Modified by:   Ming
# @Last Modified time: 2022/7/11 13:55
import logging
import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-i', '--input',
              required=True,
              type=click.Path(),
              help="The merged alleles file")
@click.option('--st',
              required=True,
              type=click.Path(),
              help="The st.tsv file")
@click.option('-q', '--query',
              required=True,
              type=click.Path(),
              help="The query file name")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put result file")
def main(input, st, query, out):
    """
    chewBBACA结果统计脚本
    """
    logging.info(f"Read the alleles info")
    info = {}
    with open(input, 'r') as IN:
        header = next(IN).strip().split('\t')
        info['total_local_num'] = len(header) - 1
        for line in IN:
            arr = line.strip().split('\t')
            if arr[0] == query:
                called_local_num = 0
                for allele in arr[1:]:
                    if allele != "N":
                        called_local_num += 1
            info["called_local_num"] = called_local_num

    logging.info(f"Read the st info")
    with open(st, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split('\t')
            res = arr[1]

    with open(out, 'w') as OUT:
        print(*["样本名称", "总位点数目", "鉴定到的位点数目", "分型"], sep="\t", file=OUT)
        print(*["Query", info['total_local_num'], info["called_local_num"], res], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
