#!/usr/bin/env python
# @CreateTime       : 2022/09/01
# @Author           : mengxf
# @version          : v1.0.1
# @LastModified     : 2022/09/01
# @ModifiedDetail   : 无
# @description      : 拆分swissprot数据库uniprot_sprot.fasta.gz到病毒、细菌、真菌字库，然后合并成宏基因组库。
#                     就是删除了真核动植物等等。
import click
import logging
from lib.make_db import MakeDB
from lib.analysis import Analysis


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

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
@click.option("-o", "--outable", default="./swissprot_result.tsv", show_default=True, help="输出结果表格.")
@click.option("--softw_align", default="blastp", show_default=True, help="选择比对软件.", 
            type=click.Choice(["diamond", "blastp"], case_sensitive=False))
@click.option("--db_select", default="Meta", show_default=True, help="数据库选择.", 
            type=click.Choice(["Meta", "Archaea", "Bacteria", "Viruses", "Fungi"]))
@click.option("-t", "--threads", default="4", show_default=True, help="线程数.")
def analysing(infa, softw_align, db_select, db_path, outable, threads):
    """swissprot分析"""
    analysing = Analysis(infa, softw_align, db_select, db_path, outable, threads)
    analysing.execute()


if __name__ == "__main__":
    cli()
