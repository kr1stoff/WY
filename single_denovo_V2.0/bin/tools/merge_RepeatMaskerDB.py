#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-11-24 15:44 
# @fuction: 合并Dfam 和 RepBase 的数据库注释文件
import os
import sys
import logging
import click
from pathlib import Path
import re
import pandas as pd
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def parse_embl(file):
    embl_dic = {}
    with open(file,'r',encoding='utf-8')as f:
        for line in f:
            if not line: continue
            if line.endswith("*"): continue
            if line.startswith("XX"): continue

            if re.search('^ID',line):
                id,length = re.findall(r"^ID   (.*?); .*; (\d+) BP\.$",line)[0][0:]
                embl_dic.setdefault(id,{})
                embl_dic[id]['length'] = int(length)

            if re.search('^DE',line):   # RepbaseID
                RepbaseID = re.findall(r"RepbaseID: (.*?)$",line)[0]
                embl_dic[id]['RepbaseID'] = RepbaseID if RepbaseID else '-'
            
            if re.search('^CC',line):
                if re.search(' Type: ',line):
                    Type = re.findall(r"Type: (.*?)$",line)[0]
                    embl_dic[id]['Type'] = Type.replace('?','') if Type else '-'
                
                elif re.search('SubType: ',line):
                    SubType = re.findall(r"SubType: (.*?)$",line)[0]
                    embl_dic[id]['SubType'] = SubType.replace('?','') if SubType else '-'
                
                elif re.search('Species: ',line):
                    Species = re.findall(r"Species: (.*?)$",line)[0]
                    embl_dic[id]['Species'] = Species if Species else '-'
                
                elif re.search('SearchStages: ',line):
                    SearchStages = re.findall(r"SearchStages: (.*?)$",line)[0]
                    embl_dic[id]['SearchStages'] = SearchStages if SearchStages else '-'
                
                elif re.search('BufferStages:',line):
                    BufferStages = re.findall(r"BufferStages:(.*?)$",line)[0]
                    embl_dic[id]['BufferStages'] = BufferStages if BufferStages else '-'
                
                elif re.search('Source: ',line):
                    Source = re.findall(r"Source: (.*?)$",line)[0]
                    embl_dic[id]['Source'] = Source if Source else '-'
        # 输出
        df = pd.DataFrame.from_dict(embl_dic, orient='index')
        df = df.reset_index().rename(columns={'index': 'ID'})
        return df


def parse_lib(file):
    lib_dic = {}
    with open(file,'r',encoding='utf-8')as f:
        for line in f:
            if not line: continue
            id,substr,Species = line.strip().split('\t')[0:3]
            lib_dic.setdefault(id,{})
            Type,subType = '-','-' 
            if re.search(r'\/',substr):
                Type,subType = re.split(r'\/',substr)[0:]
            else:
                Type = substr
            lib_dic[id]['Type'] = Type
            lib_dic[id]['SubType'] = subType
            lib_dic[id]['Species'] = Species
    # 输出
    df = pd.DataFrame.from_dict(lib_dic, orient='index')
    df = df.reset_index().rename(columns={'index': 'ID'})
    return df



@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-1','--embl', required=True,type=click.Path(exists=True),help="输入RMRB.embl")
@click.option('-2','--lib', required=True,type=click.Path(exists=True),help="输RepeatMasker.lib.txt")
@click.option('-o','--out', required=True,default='./result.txt',type=click.Path(),help="输出文件,默认[./]")
def main(embl,lib,out):
    """合并Dfam 和 RepBase 的数据库注释文件"""
    df1 = parse_embl(embl)
    df2 = parse_lib(lib)
    
    df = pd.concat([df1,df2],axis=0,join='outer')
    df.fillna('-',inplace=True)
    df.drop_duplicates(subset='ID', keep='first',inplace=True)
    df.to_csv(out,sep='\t',encoding='utf-8',index=False)
    logging.info(f"输出 {out} ")




if __name__ == '__main__':
    main()