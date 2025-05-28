#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/8/7 11:02
# @Last Modified by:   Ming
# @Last Modified time: 2023/8/7 11:02
import logging
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
sys.path.append(str(d_lib))
from Primer3 import Setting

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def merge_result(d_in: Path, info_region, f_out):
    """
    合并每个区域的设计结果

    :param d_in: The input dir of sequence
    :param info_region: The dict contain flag to region name info
    :param f_out: The output file name
    """
    res = []
    new_header = None
    for flag, name in info_region.items():
        f_in = d_in.joinpath(f"{flag}/primer.xls")
        with open(f_in, 'r') as IN:
            header = next(IN).strip().split("\t")
            new_header = header + ["TargetRegion"]
            for line in IN:
                arr = line.strip().split("\t")
                arr.append(name)
                res.append(arr)

    flag = 0
    with open(f_out, 'w') as OUT:
        print(*new_header, sep="\t", file=OUT)
        for i in res:
            i[0] = flag
            print(*i, sep="\t", file=OUT)
            flag += 1


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-r', '--region',
              required=True,
              type=click.Path(),
              help="The region file for amplified region")
@click.option('-f', '--fasta',
              required=True,
              type=click.Path(),
              help="The fasta file for design primer")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
@click.option("--config",
              required=True,
              type=click.Path(),
              help="The YAML config file for the sequence primer design")
@click.option('--primer3',
              type=click.Path(),
              help="The Primer3 path")
def main(region, fasta, out, config, primer3):
    """
    切分区域的bed文件，生成primer3所需的配置文件
    """
    f_region = Path(region).absolute()
    f_fasta = Path(fasta).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)
    f_config = Path(config).absolute()

    # Primer3 参数
    logger.info(f"Parse the config file {f_config}")
    p3_config = yaml.safe_load(open(f_config, 'r').read())

    logger.info(f"Parse the fasta file {f_fasta}")
    info_seq = {}
    for record in SeqIO.parse(f_fasta, 'fasta'):
        info_seq[record.name] = record.seq

    logger.info(f"Parse the region file: {f_region} and run the Primer3")
    flag = 0
    info_region = {}
    with open(f_region, 'r') as IN:
        for line in IN:
            arr = line.strip().split("\t")
            start, end = [int(i) for i in arr[1:3]]
            info_region[flag] = arr[3]
            d_split = d_out.joinpath(f"{flag}")
            d_split.mkdir(exist_ok=True, parents=True)
            f_setting = d_split.joinpath("setting.p3")
            with open(f_setting, 'w') as OUT:
                setting = Setting(arr[0], info_seq[arr[0]])
                # 固定参数
                for i, j in p3_config["design"].items():
                    setting.set(i, j)
                # 区域参数
                ## 排除的区域
                setting.set("SEQUENCE_EXCLUDED_REGION", f"{start + 1},{end - start}")
                ## 上下游引物区域
                pos_left = max(start - 100 + 1, 1)
                len_left = start - pos_left + 1
                pos_right = end + 2
                len_right = min(len(info_seq[arr[0]]) - end + 1, 100)
                setting.set("SEQUENCE_PRIMER_PAIR_OK_REGION_LIST", f"{pos_left},{len_left},{pos_right},{len_right}")

                print(setting, file=OUT)

            logger.debug(f"Start to run Primer3")
            f_result = d_split.joinpath("primers.txt")
            if not primer3:
                f_config = d_config.joinpath("software.yml")
                software = yaml.safe_load(open(f_config, 'r'))
                primer3 = software["primer3"]
            cmd = f"{primer3} {f_setting} > {f_result}"
            try:
                (p3out, p3error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
            except:
                sys.exit(logger.error(p3error))

            # 结果解析
            python = software["python"]
            cmd = f"{python} {d_bin}/p3_parser.py -f {f_result} -o {d_split} -p external"
            try:
                (parse_out, parse_error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
            except:
                sys.exit(logger.error(parse_error))

            flag += 1

    f_out = d_out.joinpath("primer.xls")
    logger.info(f"Merge the result to {f_out}")
    merge_result(d_out, info_region, f_out)


if __name__ == "__main__":
    main()
