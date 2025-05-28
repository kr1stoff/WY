#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/7/26

import logging
import click
import pandas as pd


#运行日志格式设置：时间+运行提示日志
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
help = dict(help_option_names=['-h', '--help'])

def parse_table(table):
    df = pd.DataFrame(pd.read_csv(table,sep='\t',encoding='utf-8',names=['taxid','kingdom','phylum','class','order','family','genus','species']))
    





#参数
@click.command(context_settings=help)
@click.option('-i','--table',required=True,type=click.STRING,help="karnken.out.lineage")
@click.option('-o', '--out',required=True,type=click.Path(),help="out file")



def main(table,out):
    """
    Parse karnken.out.lineage
    """



#调用全局函数
if __name__ == "__main__":
    main()