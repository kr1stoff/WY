#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/1/24 18:59
# @Last Modified by:   Ming
# @Last Modified time: 2022/1/24 18:59
import logging
import click
import os

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-f', '--feature',
              required=True,
              type=click.Path(),
              help="The feature table of OTU/ASV")
@click.option('-s', '--sequence',
              type=click.Path(),
              help="The represent sequence")
@click.option('-n', '--name',
              default="ASV",
              show_default=True,
              type=click.Choice(["ASV", "OTU"]),
              help="The name prefix for new name")
@click.option('-w', '--width',
              default=6,
              show_default=True,
              type=int,
              help="The number width for OTU/ASV name")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put dir")
def main(feature, sequence, name, width, out):
    """
    Rename the filter table and represent sequence
    
    *Sequence not in the feature table will abandon*
    """
    out = os.path.abspath(out)
    os.makedirs(out, exist_ok=True)
    logging.info(r"Rename the feature table")
    f_feature = os.path.join(out, "feature-table.tsv")
    info = {}
    i = 1
    with open(feature, 'r') as IN, open(f_feature, 'w') as OUT:
        for line in IN:
            if line.startswith("#"):
                print(line.strip(), file=OUT)
            else:
                arr = line.strip().split('\t')
                num = str(i).zfill(width)
                n_name = f"{name}{num}"
                info[arr[0]] = n_name
                arr[0] = n_name
                print(*arr, sep="\t", file=OUT)
                i += 1
    if sequence:
        logging.info(r"Rename the sequence file")
        f_sequence = os.path.join(out, "dna-sequences.fasta")
        flag = 0
        with open(sequence, 'r') as IN, open(f_sequence, 'w') as OUT:
            for line in IN:
                if line.startswith(">"):
                    name = line.strip()[1:]
                    if name in info:
                        flag = 1
                        print(f">{info[name]}", file=OUT)
                    else:
                        flag = 0
                else:
                    if flag:
                        print(line.strip(), file=OUT)


if __name__ == "__main__":
    main()
