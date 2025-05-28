#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/1/29 15:38
# @Last Modified by:   Ming
# @Last Modified time: 2022/1/29 15:38
import logging

import click
import pandas as pd

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-f', '--feature',
              required=True,
              type=click.Path(),
              help="The feature table")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put file")
def main(feature, out):
    """
    Get the min depth of feature table
    """
    logging.info(f"Parse the feature table")
    df = pd.read_csv(feature, sep="\t", index_col=0)
    res = int(min(df.sum()))
    with open(out, 'w') as OUT:
        print(res, file=OUT)


if __name__ == "__main__":
    main()
