#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/20 20:12
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/20 20:12
import logging
import sys
import urllib
from pathlib import Path

import click
from Bio import Entrez
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function

#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-t', '--target',
              required=True,
              type=click.Path(),
              help="The Target name list")
@click.option('--source',
              required=True,
              default="nucleotide",
              show_default=True,
              type=click.Choice(["nucleotide"]),
              help="The source database name")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put file name")
def main(target, source, out):
    """
    获取ncbi下载的数据的metadata

    TODO: 暂时只尝试过解析从nucleotide的病毒
    """
    f_target = Path(target).absolute()
    f_out = Path(out).absolute()

    logger.info(f"Parse the target {f_target}")
    ids = set()
    with open(f_target, 'r') as IN:
        for line in IN:
            ids.add(line.strip())
    ids = list(ids)

    Entrez.email = "616896512@qq.com"
    try:
        handle = Entrez.efetch(db=source, id=ids, rettype="gb", retmode="text")
    except urllib.error.HTTPError:
        sys.exit(logging.error("Can't get info through Entrez, please check your input and database"))

    res = {i: {"country": None, "collection_date": None} for i in ids}
    for chromosome in SeqIO.parse(handle, "genbank"):
        # 带有版本信息
        name = chromosome.id
        for record in chromosome.features:
            if "country" in record.qualifiers:
                res[name]["country"] = record.qualifiers["country"][0].strip().split(':')[0]
            if "collection_date" in record.qualifiers:
                res[name]["collection_date"] = record.qualifiers["collection_date"][0].strip().split('-')[-1]

    logger.info(f"Out put the result to {f_out}")
    with open(f_out, 'w') as OUT:
        print(*["Name", "Country", "Data"], sep="\t", file=OUT)
        for i, j in res.items():
            print(*[i, j["country"], j["collection_date"]], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
