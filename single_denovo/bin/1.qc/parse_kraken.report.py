#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/10/26 16:01

import click
import pandas as pd
import os 

def parse_taxid_file(file):
    taxid_list = []
    with open(file,'r',encoding='utf-8') as f:
        for line in f:
            if not line: continue
            taxid_list.append(int(line.strip()))
    return taxid_list


def parse_report(infile,taxid_list:list):
    out = os.path.dirname(os.path.abspath(infile))
    print(taxid_list)
    data = pd.read_csv(infile,
                sep='\t',
                encoding='utf-8',
                usecols=[0,1,3,4,5],
                names=['丰度占比(%)','序列数','分类级别','Taxid','物种名称'],
                dtype={"序列数": int ,"丰度占比":float})

    df = pd.DataFrame(data)
    df = df[df['分类级别'] == 'S']                   # 筛选种
    df = df[~df['Taxid'].isin(taxid_list)]          # 筛选掉taxid list 中的物种
    df['物种名称'] = df['物种名称'].map(str.strip)   # 去空格
    df.sort_values(by=["序列数"],ascending=False,inplace=True)

    new_df = df[['Taxid','物种名称','序列数','丰度占比(%)']]
    new_df.to_csv(f"{out}/abundance.txt",sep='\t',index=False,encoding='utf-8')



## Options
@click.command(context_settings = dict(help_option_names=['-h', '--help']))
@click.option('-i','--infile',required=True,type=click.Path(),help="kraken report路径")
@click.option('-l','--taxid_file',type=click.Path(),help="需要过滤掉的taxid列表")

def main(infile,taxid_file):
    """解析kraken report，生成物种丰度统计表"""
    taxid_list = parse_taxid_file(taxid_file) if taxid_file else []
    parse_report(infile,taxid_list)




if __name__ == "__main__":
    main()