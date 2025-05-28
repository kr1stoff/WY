#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/3/1 13:11
# @Last Modified by:   Ming
# @Last Modified time: 2023/3/1 13:11
import logging

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')
#### Some Function
dict_degenerate = {('A',): 'A',
                   ('C',): 'C',
                   ('G',): 'G',
                   ('T',): 'T',
                   ('A', 'G'): 'R',
                   ('C', 'T'): 'Y',
                   ('C', 'G'): 'S',
                   ('A', 'T'): 'W',
                   ('G', 'T'): 'K',
                   ('A', 'C'): 'M',
                   ('C', 'G', 'T'): 'B',
                   ('A', 'G', 'T'): 'D',
                   ('A', 'C', 'T'): 'H',
                   ('A', 'C', 'G'): 'V',
                   ('A', 'C', 'G', 'T'): 'N'}


def get_degenerate(l_alpha):
    return dict_degenerate[tuple(sorted(l_alpha))]


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-i', '--fin',
              required=True,
              type=click.Path(),
              help="The file generate by align_to_reference.py")
@click.option('-o', '--fout',
              required=True,
              type=click.Path(),
              help="The output file ")
def main(fin, fout):
    """
    DESC
    """

    with open(fin, 'r') as IN:
        ref = next(IN).strip()
        a = next(IN).strip()
        t = next(IN).strip()
        g = next(IN).strip()
        c = next(IN).strip()

    with open(fout, 'w') as OUT:
        print(ref, file=OUT)
        for i in range(len(ref)):
            tmp = set([ref[i], a[i], t[i], g[i], c[i]])
            tmp.remove('-')
            print(tmp)
            alpha = get_degenerate(list(tmp))
            print(alpha, end="", file=OUT)

    pass


if __name__ == "__main__":
    main()
