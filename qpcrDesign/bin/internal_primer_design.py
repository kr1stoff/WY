#!/usr/bin/env python

import click
from pathlib import Path
import yaml
from subprocess import run
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


@click.command()
@click.option('-f', '--fasta', type=click.Path(exists=True), help='输入用于设计引物的FASTA序列文件.')
@click.option('-o', '--outdir', help='结果输出目录.')
@click.option('-r', '--size_range', default='70-110', show_default=True, help='设计产物长度范围.')
def main(fasta, outdir, size_range):
    home = Path(__file__).parents[1]
    tmpltp3 = home.joinpath('etc/template.p3')
    primer3 = yaml.safe_load(open(home.joinpath('conf/software.yml')).read())['primer3']
    outdir = Path(outdir)
    #读FASTA
    nucbs = set(list('ATGCNRYMKSWHBVD'))
    with open(fasta) as f:
        header = next(f).strip().replace('>', '')
        bases = ''.join([line.strip() for line in f])
        if set(list(bases)).difference(nucbs): #确定是核酸序列
            raise Exception(f'输入FASTA文件异常, 请查看! {fasta}')
    #生成primer3core输入文件
    stp3 = outdir.joinpath('setting.p3')
    with open(tmpltp3) as f, open(stp3, 'w') as g:
        g.write(f'SEQUENCE_ID={header}\nSEQUENCE_TEMPLATE={bases}\nPRIMER_PRODUCT_SIZE_RANGE={size_range}\n')
        g.write(f.read())
    #命令行
    cmd = f'{primer3} < {stp3} > {outdir}/primers.txt'
    logging.info(cmd)
    run(cmd, shell=True, executable='/bin/bash')

if __name__ == '__main__':
    main()
