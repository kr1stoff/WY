#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/2/13 14:01
# @Last Modified by:   Ming
# @Last Modified time: 2023/2/13 14:01
import logging
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--internal',
              required=True,
              type=click.Path(),
              help="The internal info file")
@click.option('--external',
              required=True,
              type=click.Path(),
              help="The external info file")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The output file name")
def main(internal, external, out):
    """
    合并内外引物设计结果
    """
    f_internal = Path(internal).absolute()
    f_external = Path(external).absolute()
    f_out = Path(out).absolute()

    logger.info(f"Parse the internal file {f_internal}")
    info_internal = {}
    with open(f_internal, 'r') as IN:
        header_internal = next(IN).strip().split('\t')
        header_internal = [f"Internal_{i}" for i in header_internal]
        for line in IN:
            arr = line.strip().split("\t")
            name = f"{arr[1]}_{arr[0]}"
            info_internal[name] = arr

    logger.info(f"Parse the external file {f_external}")
    flag = 1
    with open(f_external, 'r') as IN, open(f_out, 'w') as OUT:
        header_external = next(IN).strip().split("\t")
        n_header = ["Index"] + header_internal + [f"External_{i}" for i in header_external]
        print(*n_header, sep="\t", file=OUT)
        for line in IN:
            arr = line.strip().split("\t")
            name = arr[2]
            if name in info_internal:
                print(*[flag, *info_internal[name], *arr], sep="\t", file=OUT)
                flag += 1


if __name__ == "__main__":
    main()
