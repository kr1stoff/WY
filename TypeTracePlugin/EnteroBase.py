#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/12/19 16:46
# @Last Modified by:   Ming
# @Last Modified time: 2023/12/19 16:46
import gzip
import logging
import subprocess
import sys
from pathlib import Path

import click
import yaml
from Bio import SeqIO

from lib.api import EnteroBase

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.1.0")
@click.option("--debug",
              is_flag=True,
              show_default=True,
              help="Whether turn on the debug mode")
@click.pass_context
def cli(ctx, debug):
    """
    EnteroBase工具包
    """
    if debug:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ctx.ensure_object(dict)
    ctx.obj['DEBUG'] = debug


@click.command()
@click.option("-o", "--out",
              type=click.Path(),
              help="Out put the available scheme to file")
def list(out):
    """
    列出数据库中所有可用的schemes并打印，如果指定out则将可用schemes输出到out文件
    """
    res = EnteroBase().list()
    logger.info(f"Fetch the available from EnteroBase")
    if out:
        f_out = Path(out).absolute()
        with open(f_out, 'w') as OUT:
            print(*res, sep="\n", file=OUT)
    else:
        print(*res, sep="\n")


@click.command()
@click.option("--scheme",
              required=True,
              help="The scheme name")
@click.option("-o", "--out",
              type=click.Path(),
              default="./",
              show_default=True,
              help="The output dir")
@click.option("-t", "--threads",
              default=10,
              type=int,
              show_default=True,
              help="The max thread number to use when download the allele sequence")
@click.option("--force",
              is_flag=True,
              default=False,
              show_default=True,
              help="Whether redownload the exists file")
@click.option("--fail",
              is_flag=True,
              default=False,
              show_default=True,
              help="Only download the allele sequence in the fail.txt")
def download(scheme, out, threads, force, fail):
    """
    下载目标scheme在EnteroBase的所有信息
    """
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)
    app = EnteroBase()
    if fail:
        f_fail = d_out.joinpath("fail.txt")
        logger.info(f"Start to download file in file {f_fail}")
        app.download_fail(f_fail, threads=threads)
    else:
        logger.info(f"Start to download the {scheme} to {d_out}")
        app.download(scheme, d_out, threads, force)
        logger.info(f"Job finish, please check the result before build a database")


@click.command()
@click.pass_context
@click.option("-s", "--scheme",
              type=click.Path(),
              required=True,
              help="The scheme table")
@click.option("-i", "--dir_in",
              type=click.Path(),
              required=True,
              help="The input dir with allele sequence")
@click.option("-o", "--dir_out",
              type=click.Path(),
              required=True,
              help="The output dir of the database")
@click.option("-p", "--prefix",
              required=True,
              help="The output prefix for the scheme")
def build(ctx, scheme, dir_in, dir_out, prefix):
    """
    构建EnteroBase来源的cgMLST数据库
    """
    d_in = Path(dir_in).absolute()
    d_out = Path(dir_out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)
    f_fasta = d_out.joinpath(f"{prefix}.fasta")
    f_length = d_out.joinpath(f"{prefix}.length")

    # scheme 处理
    f_scheme = d_out.joinpath(f"{prefix}.txt")
    logger.info(f"Deal the scheme info {f_scheme}")
    with gzip.open(scheme, 'rt') as IN, open(f_scheme, 'w') as OUT:
        header = next(IN).strip().split('\t')
        locis = header[1:]
        print(*header, sep="\t", file=OUT)
        for line in IN:
            print(line.strip(), file=OUT)

    # 序列处理
    logger.info(f"Merge the allele sequence to {f_fasta} and output the length info to {f_length}")
    with open(f_fasta, 'w') as OUT1, open(f_length, 'w') as OUT2:
        for loci in locis:
            f_allele = d_in.joinpath(f"{loci}.fasta.gz")
            if f_allele.exists():
                for record in SeqIO.parse(gzip.open(f_allele, 'rt'), "fasta"):
                    print(f">{record.name}\n{record.seq}", file=OUT1)
                    print(*[record.name, len(record.seq)], sep="\t", file=OUT2)
            else:
                logger.error(f"{f_allele} not exist")
                sys.exit(1)

    logger.info(f"Build the blast database")
    f_config = Path(__file__).parent.joinpath("config/software.yml").absolute()
    cfg = yaml.safe_load(open(f_config, 'r').read())
    cmd = f"{cfg['makeblastdb']} -in {f_fasta} -out {d_out}/{prefix} -dbtype nucl -input_type fasta -title '{prefix}'"
    if ctx.obj['DEBUG']:
        logger.debug(cmd)
    try:
        subprocess.run(cmd, shell=True, check=True)
    except Exception as e:
        logger.error(e)
        sys.exit(1)


cli.add_command(list)
cli.add_command(download)
cli.add_command(build)

if __name__ == "__main__":
    cli()
