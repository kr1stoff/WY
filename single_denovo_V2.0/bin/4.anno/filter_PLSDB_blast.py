#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-11-29 14:52 
# @fuction: 
import os
import sys
import logging
import click
import pandas as pd 
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def parse_tab(file):
    df = pd.read_csv(file,sep='\t',encoding='utf-8',header=0)
    df.fillna('-',inplace=True)
    return df

def parse_blast(file):
    df = pd.read_csv(file,sep='\t',encoding='utf-8',header=None)
    columns = ['qid','sid','identity','length' ,'mismatch' ,'gapopen' ,
               'qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore' ,
               'qlen' ,'slen' ,'qcovs' ,'qcovhsp' ,'qcovus']
    df.columns = columns
    filter_df = df.groupby('sid').head(1)  # 比对到同一个参考只保留1条记录
    filter_df = filter_df.drop_duplicates(subset=['qid','qstart', 'qend'], keep='first')
    return filter_df



@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-1','--tab', required=True,type=click.Path(exists=True),help="PLSDB质粒数据库注释表格")
@click.option('-2','--infile', required=True,type=click.Path(exists=True),help="blatsn结果")
def main(infile,tab):
    """blast过滤，合并注释信息"""
    dir = os.path.dirname(os.path.abspath(infile))
    tab_df = parse_tab(tab)
    blast_df = parse_blast(infile)

    df = pd.merge(blast_df,tab_df,how='left',left_on='sid',right_on='NUCCORE_ACC')
    df = df[df['NUCCORE_Completeness']=='complete']
    
    df = df.drop(columns=['NUCCORE_ACC','NUCCORE_Completeness'])
    out = os.path.join(dir,'plasmid.txt')
    df.to_csv(out,sep='\t',encoding='utf-8',index=False)
    logging.info(f'[Output]: {out}')






if __name__ == '__main__':
    main()
    