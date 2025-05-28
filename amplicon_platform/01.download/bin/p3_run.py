#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/15 10:55
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/15 10:55
import logging
import os
import subprocess
import sys
from pathlib import Path

import click
import yaml
from Bio import SeqIO

# 添加搜索路径
d_bin = Path(__file__).parent
d_base = d_bin.parent
d_lib = d_base.joinpath("lib")
d_config = d_base.joinpath("config")
sys.path.append(d_lib)
from Primer3 import Setting

# 全局设置
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--ptype',
              required=True,
              type=click.Choice(["internal","external"]),
              help="The primer type")
@click.option('-f', '--fasta',
              required=True,
              type=click.Path(),
              help="The fasta file for design primer")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put dir")
@click.option('--product_range',
              default=["120-500"],
              show_default=True,
              multiple=True,
              help="The product length range")
@click.option('--primer_opt_size',
              default=23,
              show_default=True,
              help="Optimum length of a primer")
@click.option('--primer_length_range',
              default="18-30",
              show_default=True,
              help="The primer length range")
@click.option('--primer_num',
              default=200,
              show_default=True,
              help="The primer num return")
@click.option('--primer_opt_tm',
              default=60.0,
              show_default=True,
              help="Optimum melting temperature for a primer")
@click.option('--primer_tm_range',
              default="58.0-62.0",
              show_default=True,
              help="The primer melting temperature range")
@click.option('--primer_opt_gc_percent',
              default=50,
              show_default=True,
              help="Optimum GC percent")
@click.option('--primer_gc_percent_range',
              default="40.0-60.0",
              show_default=True,
              help="The primer gc percent range")
@click.option("--config",
              help="The YAML config file for the mPCR/pPCR")
@click.option('--primer3',
              type=click.Path(),
              help="The Primer3 path")
def main(fasta, out, product_range, primer_opt_size, primer_length_range, primer_num, primer_opt_tm, primer_tm_range,
         primer_opt_gc_percent, primer_gc_percent_range, config, primer3):
    """
    Decorators for Primer3 run

    *配置文件内的优先级高于命令行参数*
    """
    f_fasta = Path(fasta).absolute()
    f_config = Path(config).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True)

    # Primer3 参数
    PRIMER_PRODUCT_SIZE_RANGE = ' '.join(product_range)
    PRIMER_OPT_SIZE = primer_opt_size
    PRIMER_MIN_SIZE, PRIMER_MAX_SIZE = primer_length_range.strip().split('-')
    PRIMER_NUM_RETURN = primer_num
    PRIMER_OPT_TM = primer_opt_tm
    PRIMER_MIN_TM, PRIMER_MAX_TM = primer_tm_range.strip().split('-')
    PRIMER_OPT_GC_PERCENT = primer_opt_gc_percent
    PRIMER_MIN_GC, PRIMER_MAX_GC = primer_gc_percent_range.strip().split('-')

    if config:
        logging.info(f"Parse the config file {f_config}")
        p3_config = yaml.safe_load(open(f_config, 'r').read())

    logging.info(f"Parse the fasta file {f_fasta}")
    f_setting = d_out.joinpath("setting.p3")
    with open(f_setting, 'w') as OUT:
        for record in SeqIO.parse(fasta, "fasta"):
            setting = Setting(record.id, record.seq)
            setting.set("PRIMER_PRODUCT_SIZE_RANGE", PRIMER_PRODUCT_SIZE_RANGE)
            setting.set("PRIMER_OPT_SIZE", PRIMER_OPT_SIZE)
            setting.set("PRIMER_MIN_SIZE", PRIMER_MIN_SIZE)
            setting.set("PRIMER_MAX_SIZE", PRIMER_MAX_SIZE)
            setting.set("PRIMER_NUM_RETURN", PRIMER_NUM_RETURN)
            setting.set("PRIMER_OPT_TM", PRIMER_OPT_TM)
            setting.set("PRIMER_MIN_TM", PRIMER_MIN_TM)
            setting.set("PRIMER_MAX_TM", PRIMER_MAX_TM)
            setting.set("PRIMER_OPT_GC_PERCENT", PRIMER_OPT_GC_PERCENT)
            setting.set("PRIMER_MIN_GC", PRIMER_MIN_GC)
            setting.set("PRIMER_MAX_GC", PRIMER_MAX_GC)
            if config:
                for i, j in p3_config.items():
                    if i.startswith("PRIMER"):
                        setting.set(i, j)
            print(setting, file=OUT)

    logging.info(f"Start to run Primer3")
    f_result = d_out.joinpath("primers.txt")
    if not primer3:
        f_config = d_config.joinpath("software.yml")
        software = yaml.safe_load(open(f_config, 'r'))
        primer3 = software["primer3"]
    cmd = f"{primer3} {f_setting} > {f_result}"
    try:
        (p3out, p3error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
    except:
        sys.exit(logging.error(p3error))


if __name__ == "__main__":
    main()
