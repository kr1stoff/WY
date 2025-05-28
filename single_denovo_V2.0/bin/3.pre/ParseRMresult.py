#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-11-27 16:21 
# @fuction: 
import os
import sys
import logging
import click
from pathlib import Path
import pandas as pd
import re
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def extract_number(s):
    match = re.search(r'\d+', s)
    return int(match.group()) if match else 0


def parse_file(file):
    res = []
    with open(file,'r',encoding='utf-8') as f:
        for line in f:
            if not line or line.startswith('#'): continue
            args = line.strip().split("\t")
            queryid,start,end,strand,kind = args[0],args[3],args[4],args[6],args[8]
            queryid = re.findall(r'^(.*?)_length',queryid)[0]
            kind = re.findall(r":(.*?)\"" ,kind)[0]
            res.append([queryid,start,end,strand,kind])
    df = pd.DataFrame(res,columns=['qid','qstart','qend','strand','class'])
    df.sort_values('qid', key=lambda x: x.apply(extract_number),inplace=True)
    return df
    

def parse_anno(anno):
    df = pd.read_csv(anno,sep='\t',encoding='utf-8')
    return df


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# @click.option('-1','--tbl', required=True,type=click.Path(exists=True),help="RepeatMasker.tbl文件")
@click.option('-2','--file', required=True,type=click.Path(exists=True),help="RepeatMasker.out文件")
@click.option('--anno', required=True,type=click.Path(exists=True),help="RepeatMasker注释文件")
@click.option('-o','--out', required=True,type=click.Path(),help="输出文件")
def main(file,anno,out):
    """ 预测结果与注释tab结合 """
    out = os.path.abspath(out)
    df = parse_file(file)
    anno_df = parse_anno(anno)
    merge_df = pd.merge(df,anno_df,left_on='class',right_on='ID',how='inner')
    merge_df = merge_df[['qid','qstart','qend','strand','class','Type','SubType','Species']]
    merge_df.to_csv(out,sep='\t',encoding='utf-8',index=False) 
    logging.info(f"输出 {out}")



if __name__ == '__main__':
    main()