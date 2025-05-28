#!/Bio/User/renchaobo/software/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: MingJia
# @Date:   2019-09-02 14:56:27
# @Last Modified by:   MingJia
# @Last Modified time: 2019-09-02 15:59:57
import logging
import os
import subprocess
from itertools import chain

import click
import yaml

#### Some Global variable
__version__ = '1.0.0'
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

BIN = os.path.dirname(__file__)
d_config = os.path.join(BIN, "../config")
f_software = os.path.join(d_config, "software.yml")
f_env = os.path.join(d_config, "env.yml")
software = yaml.safe_load(open(f_software, 'r').read())
env = yaml.safe_load(open(f_env, 'r').read())


#### Some Functions ####


def parse_group(file_in):
    res = {}
    with open(file_in, 'r') as IN:
        next(IN)
        for line in IN:
            arr = line.strip().split('\t')
            res.setdefault(arr[1], [])
            res[arr[1]].append(arr[0])
    return res


def is_filtered(data, compare_info, group_info):
    """
    Whether filter the line

    Keep the line if the sample number in one compare group are bigger half
    of the group's sample number
    """
    i = 0
    for compare in compare_info:
        sample_num = len(group_info[compare])
        half_num = int(sample_num / 2)
        if sum([value > 0 for value in data[i:i + sample_num]]) > half_num:
            return False
        i += sample_num
    return True


########################


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-i', '--input',
              required=True,
              type=click.Path(),
              help="The input table")
@click.option('-c', '--compare',
              required=True,
              help="The compare info(A-VS-B)")
@click.option('-g', '--group',
              required=True,
              help="The group info file")
@click.option('-l', '--level',
              type=int,
              default=7,
              show_default=True,
              help="The max levle to plot")
@click.option('-o', '--outdir',
              default='./',
              type=click.Path(),
              show_default=True,
              help="The output file")
def cli(input, compare, group, level, outdir):
    """
    LefSe analysis
    """
    group_info = parse_group(group)
    compares = compare.strip().split('-VS-')
    samples = list(chain(*[group_info[i] for i in compares]))
    lefse_input = os.path.join(outdir, f'{compare}.Lefse.xls')
    logging.info(f'Chosen samples are {samples}')
    with open(input, 'r') as IN, open(lefse_input, 'w') as OUT:
        header = next(IN).strip().split('\t')
        indexs = [header.index(i) for i in samples]
        line1 = ['Group']
        for i in compares:
            line1 += [i] * len(group_info[i])
        print(*line1, sep='\t', file=OUT)
        line2 = ['Samples'] + samples
        print(*line2, sep='\t', file=OUT)
        for line in IN:
            arr = line.strip().split('\t')
            clean_name = arr[0].strip().split(';', 1)
            # if len(clean_name) > 1:
            tmp = [float(arr[i]) for i in indexs]
            res = [arr[0]] + tmp
            print(*res, sep='\t', file=OUT)
            cmd = f"""set -e
source {software['activate']} {env['lefse']}
lefse-format_input.py {outdir}/{compare}.Lefse.xls {outdir}/{compare}.Lefse.in -c 1 -s 2 -o 1000000
run_lefse.py {outdir}/{compare}.Lefse.in {outdir}/{compare}.Lefse.res -f 0.9
lefse-plot_res.py {outdir}/{compare}.Lefse.res {outdir}/{compare}.Lefse.svg --format svg
{software['perl']} {BIN}/Extract_Diff_species.pl {outdir}/{compare}.Lefse.res {outdir}
tail -n +2 {outdir}/{compare}.Lefse.Diff.xls > {outdir}/{compare}.Lefse.Diff.xls.nohead
lefse-plot_cladogram.py {outdir}/{compare}.Lefse.Diff.xls.nohead {outdir}/{compare}.Lefse.cladogram.svg --labeled_stop_lev {level} --abrv_stop_lev {level} --min_point_size 0.5 --max_point_size 3 --radial_start_lev 1
rm {outdir}/{compare}.Lefse.Diff.xls.nohead
# mkdir -p {outdir}/{compare}
# lefse-plot_features.py {outdir}/{compare}.Lefse.in {outdir}/{compare}.Lefse.res {outdir}/{compare}/
{software['convert']} -density 300 {outdir}/{compare}.Lefse.svg {outdir}/{compare}.Lefse.png
{software['convert']} -density 300 {outdir}/{compare}.Lefse.cladogram.svg {outdir}/{compare}.Lefse.cladogram.png
"""
    logging.debug(cmd)
    subprocess.run(cmd, shell=True, executable=software["shell"])


if __name__ == "__main__":
    cli()
