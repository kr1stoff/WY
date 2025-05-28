#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/8/15 11:01
import logging
import os,re
import click
import numpy as np
import pandas as pd
from Bio import SeqIO
import sys

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)



class Filter():
    def __init__(self,p_fa,length,cov,prefix):
        self.p_fa   = p_fa
        self.minlen = length
        self.mincov = cov
        self.prefix = prefix
        self.out    = os.path.dirname(os.path.abspath(p_fa))
        self.list_contigs = []
        self.Length=[]
        self.new_fa = f"{self.out}/{self.prefix}.fasta"


    #读取fasta
    def read_fa(self):
        '''Make a list of Seq records from a fasta filepath''' 
        for record in SeqIO.parse(self.p_fa, 'fasta'):
            self.list_contigs.append(record)


    #过滤逻辑
    def filter_fa(self):
        #匹配contig.fasta的头,eg: >NODE_2897_length_56_cov_203.00
        cov_pattern   = re.compile("cov_([0-9.]+)")

        total_contigs = 0
        contigs_kept = 0
        kept_fa=[]

        for contig in self.list_contigs:
            total_contigs = total_contigs + 1

            #提取fa头中的覆盖度
            result = cov_pattern.search(contig.name)    

            #(1)未指定最小cov，则仅保留符合length要求的序列
            if not self.mincov:
                if len(contig) >= self.minlen:
                    contigs_kept = contigs_kept + 1
                    kept_fa.append(contig)

            #(2)提取的覆盖度不为空，保留符合cov和length要求的序列
            elif result:
                if float(result.group(1)) >= self.mincov and len(contig) >= self.minlen:
                    contigs_kept = contigs_kept + 1
                    kept_fa.append(contig)

            #(3)解析不到cov的，保留length符合要求的序列
            elif len(contig) >= self.minlen:
                logging.error(f'No coverage encoded in contig name {contig.name}')
                contigs_kept = contigs_kept + 1
                kept_fa.append(contig)

        #(4)过滤的序列写入文件
        SeqIO.write(kept_fa,self.new_fa,'fasta')

        logging.info("Starting fa: {0}\t,kept fa: {1}".format(total_contigs, contigs_kept))

#参数
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-fa', '--fasta',required=True,type=click.Path(),help="fasta路径")
@click.option('-p', '--prefix',required=True,type=click.STRING,help="Prefix of fasta,recommend use sample id ")
@click.option('-l', '--min_length',required=False,default=500,show_default=True,type=int,help="序列保留的最小长度")
@click.option('-c', '--min_cov',required=False,type=float,help="根据fasta的header标注的覆盖度过滤")


def main(fasta,min_length,min_cov,prefix):
    """
    通过序列长度，过滤组装序列\n
    eg: fasta的header例子 '>NODE_2897_length_56_cov_203.000000 \n'
    """
    project=Filter(fasta,min_length,min_cov,prefix)
    project.read_fa()
    project.filter_fa()





if __name__ == "__main__":
    main()