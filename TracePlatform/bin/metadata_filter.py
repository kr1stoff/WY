#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/4/19 14:53
# @Last Modified by:   Ming
# @Last Modified time: 2023/4/19 14:53
import logging
from pathlib import Path

import click
import pandas as pd

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--metadata',
              required=True,
              type=click.Path(),
              help="The metadata file")
@click.option('--include',
              multiple=True,
              help="The feature you want to include(eg. --include Country=China)")
@click.option('--exclude',
              multiple=True,
              help="The feature you want to exclude(eg. --include Country=China)")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The output result")
def main(metadata, include, exclude, out):
    """
    将metadata文件按照feature提供的条件进行过滤
    """
    f_metadata = Path(metadata).absolute()
    f_out = Path(out).absolute()
    includes = {i.strip().split('=')[0]: i.strip().split('=')[1] for i in include}
    excludes = {i.strip().split('=')[0]: i.strip().split('=')[1] for i in exclude}

    logger.info(f"Parse the metadata file: {f_metadata}")
    df = pd.read_csv(f_metadata)

    # 保留的feature
    for i, j in includes.items():
        df = df[df[i] == j]
    # 丢弃的feature
    for i, j in excludes.items():
        df = df[df[i] != j]

    logger.info(f"Write the result to file: {f_out}")
    df.to_csv(f_out, index=False)


if __name__ == "__main__":
    main()
