#!/home/earthtest/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/8/12 16:00

import os
import click
import pandas as pd
import logging
from Bio import SeqIO
from Bio.SeqUtils import GC
import sys

#获取需要读取的文件路径
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
PATH = os.path.dirname(os.path.abspath(__file__))   #获取本脚本路径
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


def add_base(col,dic):
    id = col[0]
    pos = int(col[1])-1
    base = dic[id]['seq'][pos]
    return base


####参数######################################################################################
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i','--depth_file', required=True,type=click.Path(exists=True),help="Depth file")
@click.option('-ref', required=True,type=click.Path(exists=True),help="Reference")
@click.option('-o','--out', required=True,type=click.Path(),help="Total stat file")
@click.option('-a','--all', required=True,type=click.Path(),help="All seq stat file")


def main(depth_file,ref,out,all):
    ref_dic={}
    logging.info(f"Parse {ref}")
    for record in SeqIO.parse(ref, 'fasta'):
        #找出N碱基位点
        N_pos=[]
        for idx, letter in enumerate(record.seq):
            if letter == 'N': N_pos.append(idx)     #从0开始

        ref_dic[record.id]={}
        ref_dic[record.id]['seq']=record.seq
        ref_dic[record.id]['length']=len(record.seq)
        ref_dic[record.id]['gc']=round((record.seq.count("G")+record.seq.count("C"))/(len(record.seq)-record.seq.count("N"))*100,2)
        ref_dic[record.id]['N']=record.seq.count("N")
        ref_dic[record.id]['N_pos']=N_pos


    logging.info(f"Parse {depth_file}")
    total_df = pd.DataFrame(pd.read_csv(depth_file,sep='\t',names=['seq','pos','depth']))
    total_df['base']=total_df.apply(add_base,axis=1,dic=ref_dic)    #axis=1 按行
    ref_id = total_df['seq'].unique().tolist()

    #总序列
    df = total_df.drop(total_df[total_df['base'] == 'N'].index)  #删除含N的行计算
    dep_med =df['depth'].median()
    cov     =f"{len(df[df['depth'] > 0])/len(df)*100:.2f}"
    av_dep  =f"{df['depth'].mean():.2f}"
    dep_5x  =f"{len(df[df['depth'] >= 5])/len(df)*100:.2f}"
    dep_10x =f"{len(df[df['depth'] >= 10])/len(df)*100:.2f}"
    dep_30x =f"{len(df[df['depth'] >= 30])/len(df)*100:.2f}"
    dep_50x =f"{len(df[df['depth'] >= 50])/len(df)*100:.2f}"
    dep_100x=f"{len(df[df['depth'] >= 100])/len(df)*100:.2f}"

    logging.info(f"Output {out}")
    with open(out,'w',encoding='utf-8')as w:
        w.write(f"覆盖度(%)\t平均深度\t中位数深度\t5x覆盖(%)\t10x覆盖(%)\t30x覆盖(%)\t50x覆盖(%)\t100x覆盖(%)\n")
        w.write(f"{cov}\t{av_dep}\t{dep_med}\t{dep_5x}\t{dep_10x}\t{dep_30x}\t{dep_50x}\t{dep_100x}\n")


    #每条contig/scaffold的指标
    OUT= open(all,'w',encoding='utf-8')
    OUT.write(f"序列名称\t序列长度\tGC含量(%)\tN碱基\t覆盖度(%)\t平均深度\t中位数深度\t5x覆盖(%)\t10x覆盖(%)\t30x覆盖(%)\t50x覆盖(%)\t100x覆盖(%)\n")
    for id in ref_id:
        tmp_df = total_df[total_df['seq']== id]
        tmp_df.index = range(len(tmp_df))                 #索引从0开始
        tmp_df=tmp_df.drop(index = ref_dic[id]['N_pos'])  #删除N碱基行

        dep_med =tmp_df['depth'].median()
        cov     =f"{len(tmp_df[tmp_df['depth'] > 0])/len(tmp_df)*100:.2f}"
        av_dep  =f"{tmp_df['depth'].mean():.2f}"
        dep_5x  =f"{len(tmp_df[tmp_df['depth'] >= 5])/len(tmp_df)*100:.2f}"
        dep_10x =f"{len(tmp_df[tmp_df['depth'] >= 10])/len(tmp_df)*100:.2f}"
        dep_30x =f"{len(tmp_df[tmp_df['depth'] >= 30])/len(tmp_df)*100:.2f}"
        dep_50x =f"{len(tmp_df[tmp_df['depth'] >= 50])/len(tmp_df)*100:.2f}"
        dep_100x=f"{len(tmp_df[tmp_df['depth'] >= 100])/len(tmp_df)*100:.2f}"
        OUT.write(f"{id}\t{ref_dic[id]['length']}\t{ref_dic[id]['gc']}\t{ref_dic[id]['N']}\t{cov}\t{av_dep}\t{dep_med}\t{dep_5x}\t{dep_10x}\t{dep_30x}\t{dep_50x}\t{dep_100x}\n")
    OUT.close()
    logging.info(f"Output {all}")

if __name__ == "__main__":
    main()