#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-12-12 17:50 
# @fuction: 
import sys
import logging
import re
from pathlib import Path
import pandas as pd
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def parse_fas(file):
    msg = []
    header = ['Gene ID','Related Genes','Virulence Factors','Description','Speices']
    with open(file,'r',encoding='utf-8') as f:
        for line in f:
            if not line: continue
            # >VFG037170(gb|WP_001081754) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii 1656-2]
            pattern = r"^>(.*?)\s\((.*?)\)\s(.*?)\[(.*?)\] \[(.*?)\]$"
            if line.strip().startswith('>'):
                res = re.findall(pattern,line.strip())
                msg.append(res[0])
    df = pd.DataFrame(msg,columns=header)
    df.fillna('-',inplace=True)
    df['VFID'] = df['Description'].str.extract('.*?\((.*?)\)')
    return df


def parse_VFS_tab(file):
    df = pd.read_csv(file,sep='\t',header=0,index_col=False)
    df.fillna('-',inplace=True)
    return df


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"python3 {sys.argv[0]} <VFDB_setB_nt.fas> <VFs.tab>")
        exit(-1)
    fafile = Path(sys.argv[1]).resolve()
    tabfile = Path(sys.argv[2]).resolve()

    dir = Path(fafile).parents[0]
    outfile = Path(dir).joinpath('merge_VFDB.tab')

    df1 = parse_fas(fafile)
    df2 = parse_VFS_tab(tabfile)
    
    df = pd.merge(df1,df2,on='VFID', how='inner')
    df.drop(['Bacteria', 'Reference','VF_Name'], axis=1,inplace=True)    
    df.to_csv(outfile,sep='\t',index=False)
    logging.info(f"[Output]: {outfile}")

