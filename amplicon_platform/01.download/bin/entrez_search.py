#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/10/19 15:00
# @Last Modified by:   Ming
# @Last Modified time: 2022/10/19 15:00
import logging
import sys
import urllib
from pathlib import Path

import click
from Bio import Entrez

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option("-d", "--database",
              default="nucleotide",
              show_default=True,
              type=click.Choice(["nucleotide", "protein", "gene", "genome", "assembly"]),
              help="The database used to search")
@click.option("--idtype",
              default="gi",
              show_default=True,
              type=click.Choice(["gi", "acc"]),
              help="The idtype to return")
@click.option("--email",
              default="MingJia@qq.com",
              show_default=True,
              help="The email address to use")
@click.option("-o", "--out",
              default="result.list",
              show_default=True,
              help="The out put file")
@click.argument("keyword", nargs=-1)
def main(database, keyword, idtype, email, out):
    """
    提供Entrez的查询功能，并获取待下载序列名称
    """
    f_out = Path(out).absolute()
    query = " AND ".join(keyword)
    logging.info(f"Start to search the {database}")
    Entrez.email = email
    try:
        handle = Entrez.esearch(db=database, term=query, retmax=1)
        record = Entrez.read(handle)
        result_count = record["Count"]
        logging.info(f"Get {result_count} record for the query {query}")
        handle = Entrez.esearch(db=database, term=query, retmax=result_count, idtype=idtype)
        record = Entrez.read(handle)
    except urllib.error.HTTPError:
        logging.error("Can't get info through Entrez, please check your input and database")
        sys.exit(1)
    logging.info(f"Write the query id to {f_out}")
    with open(f_out, 'w') as OUT:
        print(*record["IdList"], sep="\n", file=OUT)


if __name__ == "__main__":
    main()
