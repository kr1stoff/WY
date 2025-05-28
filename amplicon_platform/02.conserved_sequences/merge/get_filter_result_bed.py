#!/usr/bin/python3
# @auther: stone
# @date: 2023-08-10
import sys
import pandas as pd
import os
import click
import re


def read_taxid2name(config):
    hash = {}
    with open(config, 'r') as IN:
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            taxid = arr[1]
            name = arr[0]
            hash[taxid] = name
    return hash


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_result', required=True, type=click.Path(), help="输入文件")
@click.option('-o', '--output_result', required=True, type=click.Path(), help="输出文件")
@click.option('-s', '--start_num', default=1, type=int, help="输出文件")
@click.option('-c', '--config', required=True, type=click.Path(), help="taxid和拉丁名的对应表")
@click.option('-t', '--taxid', required=True, type=click.Path(), help="taxid")
def main(input_result, output_result, config, taxid, start_num):
    hash = {}
    hash = read_taxid2name(config)
    sci_name = hash[taxid]
    with open(input_result, 'r') as IN, open(output_result, 'w') as OUT:
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            # NW_002196669.1:167340-168414_sliding:501-1000
            pattern_gene = r'(\S+):(\d+)-(\d+)_sliding:(\d+)-(\d+)'
            match_gene = re.search(pattern_gene, name)
            if match_gene:
                chr = match_gene.group(1)
                resion_start = int(match_gene.group(2))
                p1 = int(match_gene.group(4))
                p2 = int(match_gene.group(5))
                start = resion_start + p1 - 1
                end = resion_start + p2 - 1
                print(f"{chr}\t{start}\t{end}\t{sci_name}-{start_num}", file=OUT)
                start_num += 1
            # NC_001806.2_sliding:146701-146900
            pattern_genome = r'(\S+)_sliding:(\d+)-(\d+)'
            match_genome = re.search(pattern_genome, name)
            if match_genome:
                chr = match_genome.group(1)
                start = int(match_genome.group(2))
                end = int(match_genome.group(3))

                print(f"{chr}\t{start}\t{end}\t{sci_name}-{start_num}", file=OUT)
                start_num += 1


if __name__ == '__main__':
    main()
