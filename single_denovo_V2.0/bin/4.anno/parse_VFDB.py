#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-12-13 15:31 
# @fuction: 
import os
import sys
import logging
import click
from pathlib import Path
import pandas as pd
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def parse_blast_out(file):
    """
    解析blast结果
    """
    columns = ['Contig','Gene ID','identity','length' ,'mismatch' ,'gapopen' ,
               'qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore' ,
               'qlen' ,'slen' ,'qcovs' ,'qcovhsp' ,'qcovus']
    df = pd.read_csv(file,sep='\t',encoding='utf-8',header=None)
    df.columns = columns
    # df = df.groupby('Contig').head(1)  # 比对到同一个参考只保留1条记录
    df = df.drop_duplicates(subset=['Gene ID','qstart', 'qend'], keep='first')
    #条件筛选，覆盖度 
    df = df.loc[(df["qcovs"] > 80)]
    return df


def parse_tab(file):
    df = pd.read_csv(file,sep='\t',encoding='utf-8',header=0)
    df.fillna('-',inplace=True)
    return df


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-1','--blast', required=True,type=click.Path(),help="输入blast结果文件")
@click.option('-2','--depcov', required=False,type=click.Path(),help="输入深度覆盖度文件")
@click.option('--ref', required=True,type=click.Path(exists=True),help="参考序列tab")
def main(blast,depcov,ref):
    """
    解析毒力分析blast结果，并与覆盖度文件&参考注释文件合并
    """
    outdir = Path(blast).resolve().parents[0]
    outfile = Path(outdir).joinpath('Virulence.txt')

    if os.path.exists(blast):
        blast_df = parse_blast_out(blast)
        ref_df = parse_tab(ref)

        if depcov:
            depcov_df = parse_tab(depcov)
            depcov_df["Depth"] = depcov_df['Depth'].apply(lambda x: '{:.2f}'.format(x))
            depcov_df["Coverage"] = depcov_df['Coverage'].apply(lambda x: '{:.2f}'.format(x))
            total_df = pd.merge(blast_df,depcov_df,how='left',on="Contig")
            total_df = pd.merge(total_df,ref_df,how='left',on="Gene ID")
        else:
            total_df = pd.merge(blast_df,ref_df,how='left',on="Gene ID")
        
        total_df = total_df.drop_duplicates(subset=['Related Genes','Virulence Factors'], keep='first')
        total_df.fillna('-',inplace=True)
        total_df.to_csv(outfile,sep='\t',encoding='utf-8',index=False)
        logging.info(f"[Output]: {outfile}")



if __name__ == '__main__':
    main()