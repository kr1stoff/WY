#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/2/16 11:01

import os
import click
import pandas as pd
import sys



#解析文本，返回标题，和一个字典
def parse_file(file):
    with open(file, 'r') as f:
        #格式： contigs	{'marker lineage': 'k__Bacteria'
        for line in f:
            line=line.replace("# ","")
            text=line.split('\t',1)[1]           
            #eval 字符串转字典
            dic=eval(text)           
    return dic





@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--in_file',required=True,type=click.Path(),help="path of the bin_stats_ext.tsv")
@click.option('-o', '--out_file',required=True,type=click.Path(),help="path of outfile")

#################################################################################
#主函数
def main(in_file,out_file): #输入几个参数，就至少写入几个变量
    """Parse the checkm 'bin_stats_ext.tsv'"""
    
    IN = os.path.abspath(in_file)
    OUT= os.path.abspath(out_file)

    #调用函数,及输出的变量
    dic= parse_file(IN)

    #dic 转 dataframe
    df= pd.DataFrame([dic])
    total_df=df[['marker lineage','genomes','markers','Completeness','Contamination','Coding density']]
    total_df=total_df.copy()
    

    #转百分数
    total_df['Completeness']=total_df.loc[:,'Completeness'].apply(lambda x: format(x/100, '.1%')) 
    total_df['Contamination']=total_df.loc[:,'Contamination'].apply(lambda x: format(x/100, '.1%')) 
    total_df['Coding density']=total_df.loc[:,'Coding density'].apply(lambda x: format(x, '.1%')) 

    total_df.rename(columns={"marker lineage": "Marker_lineage",
                            "genomes": "参考基因组数目", 
                            "markers":"标记基因数目",
                            "Completeness": "基因组完整度",
                            "Contamination": "基因组污染度",
                            "Coding density": "cds区密度"},inplace=True)

    total_df.to_csv(OUT,sep='\t',index=False,encoding='utf-8')


               

# #调用全局函数
if __name__ == "__main__":
    main()