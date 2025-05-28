#!/usr/bin/env python
# @CreateTime       : 2022/09/01
# @Author           : mengxf
# @version          : v1.0.1
# @LastModified     : 2022/09/01
# @description      : 拆分swissprot数据库uniprot_sprot.fasta.gz到病毒、细菌、真菌子库，然后合并成宏基因组库。
#                     就是删除了真核动植物等等。
import click
import logging
from lib.make_db import MakeDB
from lib.analysis import Analysis
import sys

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


@click.group()
def cli():
    """UniProt 基因功能注释程序."""


@cli.command()
@click.option("-i", "--infa", required=True, help="输入 FASTA 数据库文件. 例:uniprot_sprot.fasta.gz.")
@click.option("-o", "--outdir", default=".", help="输出数据库目录.")
def make_db(infa, outdir):
    """根据输入的FASTA数据库文件, 拆分数据库为细菌/真菌/病毒子库和合并起来的META库."""
    makedb = MakeDB(infa, outdir)
    makedb.execute()



@cli.command()
@click.option("-i", "--infa", required=True, help="输入 FASTA Query 输入文件.")
@click.option("-d", "--db_path", required=True, help="数据库路径, make_db 创建好的库.")
@click.option("-o", "--outable", default="./swissprot_result.tsv", help="输出结果表格. [default: ./swissprot_result.tsv]")
@click.option("--softw_align", default="blastp", help="选择比对软件. [default: blastp]", type=click.Choice(["diamond", "blastp"],case_sensitive=False))
@click.option("--db_select", default="Meta", help="数据库选择. [default: Meta]", type=click.Choice(["Meta", "Archaea", "Bacteria", "Viruses", "Fungi"]))

def analysing(infa, softw_align, db_select, db_path, outable):
    """swissprot分析"""
    analysing = Analysis(infa, softw_align, db_select, db_path, outable)
    analysing.execute()


if __name__ == "__main__":
    cli()