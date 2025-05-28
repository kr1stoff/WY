#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import re
import click
import pandas as pd

def parse_eggnog(infile,desc_file,level_file):
    dir = os.path.dirname(os.path.abspath(infile))
    out = os.path.join(dir,'KEGG')
    os.makedirs(out,exist_ok=True)

    df = pd.DataFrame(pd.read_csv(infile,sep='\t',encoding='utf-8' ,comment='#',usecols=[0,11],names=['基因ID','KEGG_ko']))
    print(f"Init infile(emapper.annotations) len = {len(df)}")
    df.drop(df[df['KEGG_ko'] == '-'].index, inplace=True)    #去掉空值'-'行
    print(f"Drop empty row result,remain {len(df)}")
    df.drop_duplicates(keep='last',inplace=True)    #去重
    print(f"Drop duplicates row result,remain {len(df)}")
    df['KEGG_ko']=df['KEGG_ko'].map(lambda x:x.replace('ko:','').split(','))        #拆分行，如“ko:K16087, ko:K16089”
    df_new=df.explode('KEGG_ko')    #行，为列表才能拆分
    print(f"Split row result (included more ko) ,remain {len(df_new)}")

    # 注释文件vlookup
    desc_df = pd.DataFrame(pd.read_csv(desc_file,sep='\t',encoding='utf-8',names=['KEGG_ko','功能描述']))
    level_df = pd.DataFrame(pd.read_csv(level_file,sep='\t',encoding='utf-8',names=['KEGG_ko','一级分类','二级分类']))

    total_df = pd.merge(df_new,desc_df,how='left',on="KEGG_ko")
    total_df = pd.merge(total_df,level_df,how='left',on="KEGG_ko")
    total_df.dropna(inplace=True)
    total_df.to_csv(f"{out}/KEGG_anno.txt",sep='\t',index=False,encoding='utf-8')

    count_df = total_df.groupby(["一级分类", "二级分类"], as_index=False)['基因ID'].count()
    count_df = count_df.rename(columns={'基因ID': '基因数量'})
    count_df.to_csv(f"{out}/KEGG_anno.stat.txt",sep='\t',index=False,encoding='utf-8')


# options 
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--infile',required=True,type=click.Path(),help="emapper.annotations file")
@click.option('-db','--database',required=False,type=click.Path(),help="kegg database")


def main(infile,database):
    desc_file = os.path.join(database,'ko.txt.desc')
    level_file = os.path.join(database,'ko.txt.level')

    parse_eggnog(infile,desc_file,level_file)



if __name__ == "__main__":
    main()