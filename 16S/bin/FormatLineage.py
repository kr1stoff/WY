#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/1/26 19:24
# @Last Modified by:   Ming
# @Last Modified time: 2022/1/26 19:24
import logging
import sys

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
              help="The input file contain the lineage info")
@click.option('-n', '--name',
              required=True,
              type=click.Path(),
              help="The column name of lineage")
@click.option('--old_sep',
              default=";",
              show_default=True,
              help="The old sep for lineage")
@click.option('--new_sep',
              default="|",
              show_default=True,
              help="The new sep for lineage")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put file")
def main(input, name, old_sep, new_sep, out):
    """
    Format the lineage info in the table
    """
    with open(input, 'r') as IN, open(out, 'w') as OUT:
        header = next(IN).strip().split('\t')
        print(*header, sep="\t", file=OUT)
        try:
            index = header.index(name)
        except ValueError:
            logging.error(f"{name} not in the header of {input}")
            sys.exit(1)
        for line in IN:
            arr = line.strip().split("\t")
            name = arr[index]
            n_name = name.strip().split(old_sep)
            while n_name[-1].endswith("__"):
                if n_name[-1] == "__":
                    n_name.pop()
                else:
                    n_name[-1] += "Unclassified"
            for i in range(len(n_name)):
                if n_name[i].endswith("__"):
                    n_name[i] += "Unclassified"
            arr[index] = new_sep.join([i.strip() for i in n_name])
            print(*arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
