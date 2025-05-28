#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/7/6 14:34
# @Last Modified by:   Ming
# @Last Modified time: 2022/7/6 14:34
from Bio import SeqIO
import logging
import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('--fasta',
              required=True,
              type=click.Path(),
              help="The YAML config file for metabolome")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The out put merge file")
@click.option('-n', '--nnumber',
              default=50,
              show_default=True,
              type=int,
              help="The gap n number add to seq")
def main(fasta, out, nnumber):
    """
    如果 fasta 文件包含多个序列，则将其用 N 合并为一条，如果只有一条则仅仅重命名
    """
    logging.info(f"Parse sequence {fasta}")
    res = []
    for record in SeqIO.parse(fasta, "fasta"):
        res.append(str(record.seq))

    with open(out, 'w') as OUT:
        if len(res) > 1:
            gap = 'N' * nnumber
            new_seq = gap.join(res)
            print(f">Query\n{new_seq}", file=OUT)
        else:
            print(f">Query\n{res[0]}", file=OUT)


if __name__ == "__main__":
    main()
