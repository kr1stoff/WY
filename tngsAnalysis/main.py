#!/usr/bin/env python
# @CreateTime       : 2022/12/28
# @Author           : mengxf
# @version          : v1.0.0
# @LastModified     : 2022/12/29
# @description      : 扩增子平台分析工具集

import click
import logging
from libs.dimer_assessment import DimerAssess, DimerAssessPE
from libs.amplicon_abundance import AmpliconAbundance


#运行日志
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

@click.group(context_settings={'help_option_names':['-h','--help']})
def cli():
    """WY扩增子平台结果统计."""
@cli.command()
@click.option('-i', '--fastq1', required=True, type=click.Path(exists=True), help='输入FASTQ文件, 如果双端就是FASTQ1.')
@click.option('-I', '--fastq2', type=click.Path(exists=True), help='双端测序的FASTQ2.')
@click.option('-p', '--primer_fasta', required=True, type=click.Path(exists=True), help='输入引物FASTA文件.')
@click.option('-o', '--outdir', default='.', show_default=True, help='输出结果文件夹.')
@click.option('-P', '--prefix', default='prefix', show_default=True, help='输出结果前缀/样本名.')
def dimer_assessment(fastq1, fastq2, primer_fasta, outdir, prefix):
    """引物二聚体评估."""
    if fastq2:
        dimer = DimerAssessPE(fastq1, fastq2, primer_fasta, outdir, prefix)
    else:
        dimer = DimerAssess(fastq1, primer_fasta, outdir, prefix)
    dimer.execute()


@cli.command()
@click.option('-i', '--fastq1', required=True, type=click.Path(exists=True), help='输入下机FASTQ文件, 双端的READ1或单端.')
@click.option('-I', '--fastq2', type=click.Path(exists=True), help='输入下机FASTQ文件, 双端的READ2.')
@click.option('-r', '--reference', required=True, type=click.Path(exists=True), help='Panel参考基因组.')
@click.option('-b', '--bed', required=True, type=click.Path(exists=True), help='输入扩增子BED文件.')
@click.option('-o', '--outdir', default='.', show_default=True, help='输出结果文件夹.')
@click.option('-P', '--prefix', default='prefix', show_default=True, help='输出结果前缀/样本名.')
@click.option('-n', '--dryrun', is_flag=True, help='不执行程序,只生成shell文件.')
def amplicon_abundance(fastq1, fastq2, reference, bed, outdir, prefix, dryrun):
    """扩增子丰度分析."""
    abundance = AmpliconAbundance(fastq1, fastq2, reference, bed, outdir, prefix, dryrun)
    abundance.execute()

if __name__ == '__main__':
    cli()
