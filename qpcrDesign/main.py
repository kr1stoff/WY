#!/usr/bin/env python

# @CreateTime       : 2023/5/19
# @Author           : mengxf
# @version          : v1.3
# @LastModified     : 2023/11/15
# @version          : v1.3
# @LastModified     : 2023/11/15
# @Description      : 基于primer3的qPCR设计脚本, 替代Windows版PrimerExpress3.0.1实现自动化引物&探针设计. 
# @From             : 脚本基础来自 ChaoboRen/TracePlatform 

import logging
import sys
from pathlib import Path
import click
# custom
from lib.qPCR import Qprimer
from lib.Batch import BatchQPCR
from lib.prePCR import prePCR
from lib.degenerate import DegenerateMSA


logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.group()
def cli():
    """qPCR引物设计软件合集."""


@cli.command()
@click.option('-f', '--fasta', required=True, type=click.Path(exists=True), help='输入基因区域FASTA文件.')
@click.option('-o', '--out', default='primer3InnerDesign', show_default=True, help="引物设计结果输出目录.")
@click.option('-o', '--out', default='primer3InnerDesign', show_default=True, help="引物设计结果输出目录.")
@click.option('-r', '--size_range', default='70-110', show_default=True, help='设计产物长度范围.')
@click.option('-n', '--dryrun', is_flag=True, help='不执行程序,只生成shell文件.')
def qpcr(out, fasta, size_range, dryrun):
    """单区域qPCR设计流程, 输入GeneRegionFASTA, 输出设计结果."""
    fasta = Path(fasta).absolute()
    d_out = Path(out).absolute()
    project = Qprimer(fasta, d_out, dryrun)
    project.format()
    project.internal_design(size_range=size_range)
    project.finish()


@cli.command()
@click.option('-r', '--reference', required=True, type=click.Path(exists=True), help='输入考基因组FASTA文件.')
@click.option('-b', '--bed', required=True, type=click.Path(exists=True), help='输入保守区域BED文件.')
@click.option('-o', '--out', default='primer3InnerDesignBatch', show_default=True, help="批量设计结果输出目录.")
def qpcr_batch(reference, bed, out):
    """批量qPCR设计, 输入参考基因组FASTA和保守区域BED文件."""
    mypath = Path(__file__).resolve() #本脚本路径
    Path(out).mkdir(exist_ok=True, parents=True)
    batch = BatchQPCR(mypath, reference, bed, out)
    batch.get_region_fasta()
    batch.design1()
    batch.design2()
    batch.get_exp_txt()


@cli.command()
@click.option('-r', '--ref', type=click.Path(exists=True), help='参考基因组.')
@click.option('-b', '--bed', type=click.Path(exists=True), help='输入BED文件, 内引物评估脚本生成, 与参考基因组ACCESSION号一致.')
@click.option('-o', '--outdir', default='PrimerPlexOuterDesign', show_default=True, help='输出目录.')
@click.option('-o', '--outdir', default='PrimerPlexOuterDesign', show_default=True, help='输出目录.')
@click.option('-R', '--retry', type=click.INT, default=3, show_default=True, 
              help='如果NGS-PrimerPlex设计失败,重新筛选BED再运行的重试次数.')
def prepcr(ref, bed, outdir, retry):
    """
    巢式PCR预扩增外引物设计.
    
    BED文件格式

    NZ_CP025070.1  1365699  1365768  Brdtll prprtsss-grp_1481-1  1  B  W
    NZ_CP025070.1  1365700  1365769  Brdtll prprtsss-grp_1481-2  1  B  W
    NZ_CP025070.1  1365700  1365769  Brdtll prprtsss-grp_1481-3  1  B  W
    NZ_CP025070.1  1365700  1365769  Brdtll prprtsss-grp_1481-4  1  B  W
    """
    prepcr = prePCR(ref, bed, outdir, retry)
    prepcr.execute()


@cli.command()
@click.option('-p', '--primertable', required=True, type=click.Path(exists=True), 
              help='引物表格文件,包含内左右引物/探针/外左右引物/参考登录号.')
@click.option('-a', '--allfna', required=True, type=click.Path(exists=True), help='单物种目录下的all.fna文件.')
@click.option('-d', '--outdir', default='degenerate', show_default=True, help='输出结果目录.')
def translate_degenerate_base(primertable, allfna, outdir):
    """简并碱基转换."""
    pipe = DegenerateMSA(primertable, allfna, outdir)
    pipe.degenerate()


if __name__ == '__main__':
    cli()
