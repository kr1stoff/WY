#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-11-28 14:26 
# @fuction: 
import os
import sys
import logging
import click
from pathlib import Path
import re
import pandas as pd
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def extract_number(s):
    match = re.search(r'\d+', s)
    return int(match.group()) if match else 0

def parse_Trf_Dat(file,dir):
    df = pd.DataFrame()
    with open(file,'r',encoding='utf-8') as f:
        for line in f:
            if not line or line.startswith('#'): continue
            if line.startswith('@'):
                # @NODE_1_length_377820_cov_338.275730
                name = re.findall(r"^@(.*?)_length",line)[0]
            else:
                # TRF\tStart\tEnd\tPeriod Size\tCopy Number\tConsensus Size\tPercent Matches\tPercent Indels\tScore\tA\tC\tG\tT\tEntropy(0-2)\tleft\treight
                res = re.split(r'\s',line)[0:13]
                res = [name] + res
                df = pd.concat([df,pd.DataFrame(res).T],axis=0,join='outer')
    columns = ['ID','Start','End','Period Size','Copy Number','Consensus Size',
            'Percent Matches','Percent Indels','Score','A','G','C','T','Entropy(0-2)']
    df.columns = columns
    OUT = os.path.join(dir,'TandemRepeat.txt')  # 串联重复序列
    df.sort_values('ID', key=lambda x: x.apply(extract_number),inplace=True)
    df.to_csv(OUT,sep='\t',encoding='utf-8',index=False)
    logging.info(f"输出 {OUT}")


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i','--dat', required=True,type=click.Path(exists=True),help="trf.dat")
def main(dat):
    """解析 trf 软件dat结果"""
    dir = os.path.dirname(os.path.abspath(dat))
    parse_Trf_Dat(dat,dir)


if __name__ == '__main__':
    main()