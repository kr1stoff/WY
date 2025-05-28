#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/16 11:01
from encodings import utf_8_sig
import logging
import os,re
import click
import json
import pandas as pd
import xlwt
import sys

#运行日志格式设置：时间+运行提示日志
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
help = dict(help_option_names=['-h', '--help'])

###################################################################################################

#解析文本，返回标题，和一个列表
def parse_file(file):
    dict={}
    no_need=("contigs","bases")
    with open(file, 'r') as f:
        
        next(f)
        for line in f:
            line.strip()
            if line.startswith(no_need):
                continue
            
            li=line.strip().split(":")
            dict[li[0]]=li[1]
 
    return dict

#封装程序  main(必须参数，说明)
@click.command(context_settings=help)
@click.option('-i', '--in_file',required=True,type=click.Path(),help="path of the quast.tsv")
@click.option('-o', '--out',required=True,type=click.Path(),help="path of outfile")

#################################################################################
#主函数
def main(in_file,out): #输入几个参数，就至少写入几个变量
    #输入usage提示用途信息
    """
    Parse the quast 'predict.txt'
    """
    
    IN = os.path.abspath(in_file)

    #调用函数,及输出的变量
    logging.info(f"Parse the file {IN}")
    li= parse_file(IN)

    df= pd.DataFrame(list(li.items()))  #字典转dataframe，列标题不做index
    df.columns=["预测类型","数量"]
    
    
    df.to_csv(out,sep='\t',index=False,encoding='utf_8_sig')
    #df.to_excel(f"{out_dir}/predict_kind.xlsx",index=False,encoding='utf_8_sig')  

    logging.info(f"output the result {out}")

        

# #调用全局函数
if __name__ == "__main__":
    main()

