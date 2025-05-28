#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/14 15:37
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/14 15:37
import logging
import os
import sys

import click

BIN = os.path.dirname(__file__)
lib_path = os.path.join(BIN, "../lib")
sys.path.append(lib_path)

from Primer3 import Parser

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-f', '--file',
              required=True,
              type=click.Path(),
              help="The Primer3 out put result")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put dir")
@click.option('-p', '--ptype',
              default="internal",
              type=click.Choice(["internal", "external"]),
              help="The primer type for p3 file")
def main(file, out, ptype):
    """
    Primer3 Parser
    """
    logging.info(f"Parse the file {file}")
    res = Parser(ptype)
    logging.debug(file)
    res.load_file(file)
    f_tsv = os.path.join(out, "primer.xls")
    f_fasta = os.path.join(out, "primer.fasta")
    res.to_csv(f_tsv, sep="\t")
    res.to_fasta(f_fasta)


if __name__ == "__main__":
    main()
