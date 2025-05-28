#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-12-04 17:16 
# @fuction: 
import os
import sys
import logging
import click
from pathlib import Path
import pandas as pd
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


class ParseResfinder():
    def __init__(self,**kwargs):
        self.db_file = os.path.abspath(kwargs['db_file'])
        self.file = os.path.abspath(kwargs['infile'])
        self.flag = True if os.path.exists(self.file) else False
        self.out = os.path.dirname(self.file)

    @staticmethod
    def parse_file(file):
        df = pd.read_csv(file,sep='\t',encoding='utf-8',comment='#')
        df.fillna('-',inplace=True)
        return df

    @staticmethod
    def merge_df(df1,df2,key):
        df = pd.merge(df1, df2, on=key, how='inner')
        return df

    def run(self):
        OUT = os.path.join(self.out,'Resfinder_anno.txt')
        if self.flag:
            anno_df = self.parse_file(self.file)
            db_df = self.parse_file(self.db_file)
            db_df['Accession no.'] = db_df['Gene_accession no.'].str.extract(r'_.*?_(.*?)$')
            merge_df = self.merge_df(anno_df,db_df,key='Accession no.')
            del merge_df['Notes']
            
            logging.info(f"[Output]: {OUT}")
            merge_df.to_csv(OUT,sep='\t',encoding='utf-8',index=False)
        else:
            with open(OUT,'w',encoding='utf-8') as w:
                w.write(f"Resistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\tPosition in reference\tContig\tPosition in contig\tPhenotype_x\tAccession no.\tGene_accession no.\tClass\tPhenotype_y\tPMID\tMechanism of resistance\tRequired_gene")


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i','--infile', required=True,type=click.Path(),help="输入文件")
@click.option('-d','--db', required=True,type=click.Path(exists=True),help="数据库目录")
def main(infile,db):
    """
    解析resfinder输出文件
    """
    project = ParseResfinder(
        infile = infile,
        db_file = db
    )
    project.run()





if __name__ == '__main__':
    main()
