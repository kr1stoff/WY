#!/usr/bin/python3
# @auther: zhuzi
# @date: 2022-12-05
import re
import sys
import pandas as pd
import os
import click
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M', stream=sys.stdout)


def parse_chrome_tsv(table):
    chr_dic = {}
    with open(table, 'r', encoding='utf-8') as f:
        for line in f:
            if not line:
                continue
            lines = line.strip().split('\t')
            for chr in lines[2].strip().split(','):
                chr_dic[chr] = lines[1]
    return chr_dic


# == 计算基因数量 =====================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("blast_out", required=True, type=click.Path())
@click.option('-i', '--chr_table', required=True, type=click.Path(), help="assembly_chromosome.tsv")
@click.option('-o', '--out', required=True, type=click.Path(), help="输出文件")
def main(blast_out, out, chr_table):
    """统计比对到assembly_chromosome.tsv中的染色体的blast记录"""
    OUT = open(out, 'w', encoding='utf-8')
    chr_dic = parse_chrome_tsv(chr_table)
    with open(blast_out, 'r', encoding='utf-8') as f:
        for line in f:
            if not line:
                continue
            if line.startswith("qaccver"):
                continue
            lines = line.strip().split('\t')
            msg = '\t'.join(lines[2:])
            gene_name = chr_dic[lines[1]]
            qaccver = lines[0]
            OUT.write(f"{qaccver}\t{gene_name}\t{msg}\n")
    logging.info(f"[Output]: {out}")


if __name__ == '__main__':
    main()
