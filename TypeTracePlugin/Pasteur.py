#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/12/20 14:42
# @Last Modified by:   Ming
# @Last Modified time: 2023/12/20 14:42
import json
import logging
from pathlib import Path

import click

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("--debug",
              is_flag=True,
              show_default=True,
              help="Whether turn on the debug mode")
@click.pass_context
def cli(ctx, debug):
    """
    Pasteur数据库下载工具
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
              default="schemes.json",
              show_default=True,
              help="The out put json file")
def list(out):
    """
    列出数据库中所有可用的schemes并将其保存到输出目录的schemes.json文件中

    该数据库会有较多的scheme不在其返回的json中，可通过访问物种具体的网页获取database name，类似https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_diphtheria_seqdef
    中的数据库名称pubmlst_diphtheria_seqdef，通过构建API地址https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/schemes获取可用的schemes
    """
    from lib.api import Pasteur
    logger.info(f"Start to get available schemes from Pasteur")
    res = Pasteur().list_scheme()
    f_out = Path(out).absolute()
    f_out.parent.mkdir(exist_ok=True, parents=True)
    with open(f_out, 'w') as OUT:
        json.dump(res, OUT, indent=2)


@click.command()
@click.option("--database",
              required=True,
              help="The database name")
@click.option("--scheme_id",
              required=True,
              help="The scheme id ")
@click.option("-o", "--out",
              type=click.Path(),
              default="./",
              show_default=True,
              help="The output dir")
@click.option("-p", "--prefix",
              required=True,
              help="The output prefix for the scheme")
@click.option("-t", "--threads",
              default=10,
              type=int,
              show_default=True,
              help="The max thread number to use when download the allele sequence")
@click.option("--fail",
              is_flag=True,
              default=False,
              show_default=True,
              help="Only download the allele sequence in the fail.txt")
@click.option("--force",
              is_flag=True,
              default=False,
              show_default=True,
              help="Whether download the exists again")
def download(database, scheme_id, out, prefix, threads, fail, force):
    """
    下载目标scheme及其等位基因序列\n

    已炭疽杆菌为例：\n
    使用python main.py list下载的schemes.json文件中其信息为："B. anthracis cgMLST": "http://rest.pubmlst.org/db/pubmlst_bcereus_seqdef/schemes/2"\n
    其database name为pubmlst_bcereus_seqdef\n
    其scheme id为2\n
    下载命令为: python main.py download --database pubmlst_bcereus_seqdef --scheme_id 2\n
    """
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)
    if fail:
        from lib.api import Pasteur
        f_fail = d_out.joinpath("fail.txt")
        logger.info(f"Start to download the info from {f_fail} again")
        Pasteur().download_fail(f_fail, threads, True)
    else:
        from lib.api import Pasteur
        logger.info(f"Start to download {database} scheme {scheme_id}")
        Pasteur().download_scheme(database, scheme_id, d_out, prefix, threads, force)


cli.add_command(list)
cli.add_command(download)

if __name__ == "__main__":
    cli()
