#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/23 10:35
# @Last Modified by:   Ming
# @Last Modified time: 2022/6/23 10:35
from glob import glob
import subprocess
import logging
import click
import os

__version__ = '1.0.0'
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-i', '--input',
              required=True,
              type=click.Path(),
              help="The IRMA assemble dir")
@click.option('-o', '--output',
              required=True,
              type=click.Path(),
              help="The IRMA assemble dir")
@click.option('--db',
              required=True,
              type=click.Path(),
              help="The Flu database")
@click.option('--blastn',
              required=True,
              type=click.Path(),
              help="The blastn path")
@click.option('--param',
              required=True,
              type=click.Path(),
              help="The blastn parameter")
def main(input, output, db, blastn, param):
    """
    流感流程中IRMA结果序列比对
    """
    logging.info("Start Blast Flu all database")
    for f_fasta in glob(f"{input}/*.fasta"):
        if os.path.exists(f_fasta):
            name = os.path.basename(f_fasta).strip().split('.')[0]
            prefix = os.path.join(output, f'{name}')
            cmd = f"""{blastn} {param} -query {f_fasta} -db {db} -out {prefix}.blast.xls
cut -f 1-5 {prefix}.blast.xls | sed 1i"Query\tSubject\tIdent\tLength\tMismatch" > {prefix}.annot.xls"""
        else:
            cmd = f"""sed 1i"Query\tSubject\tIdent\tLength\tMismatch\nNA\tNA\tNA\tNA\tNA" > {prefix}.annot.xls"""
        subprocess.run(cmd, shell=True)


if __name__ == "__main__":
    main()
