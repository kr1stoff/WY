#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-12-05 15:02 
# @fuction: 
import os
import sys
import logging
import click
from pathlib import Path
import pandas as pd
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s - %(message)s', 
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


class Rgi():
    def __init__(self,**kwargs):
        self.file = os.path.abspath(kwargs['file'])
        self.depcov = os.path.abspath(kwargs['depcov'])
        self.ref = os.path.abspath(kwargs['ref'])
        
        self.out = os.path.dirname(self.file)
        self.table = os.path.join(self.out,'CARD_anno.txt')
        self.outfa = os.path.join(self.out,'CARD_contig.fa')
        self.outreffa = os.path.join(self.out,'CARD_ref.fa')


    def parse_drug_anno(self):
        self.df = pd.read_csv(self.file,sep='\t',encoding='utf-8')
        self.df.fillna('-',inplace=True)
        self.df['Contig'] = self.df['Contig'].str.strip().str.rsplit('_',n=1).str.get(0)
        
        # 输出fa
        logging.info(f"[Output]: {self.outfa}")
        with open(self.outfa,'w',encoding='utf-8') as w:
            for index, row in self.df.iterrows():
                w.write(f">{row['Contig']} {row['Start']}-{row['Stop']} ARO:{row['ARO']}|{row['AMR Gene Family']}\n")
                w.write(f"{row['Predicted_DNA']}\n")
        
        # 输出
        self.df = self.df.drop(columns = ['ORF_ID','Cut_Off','SNPs_in_Best_Hit_ARO','Best_Identities',
                        'Other_SNPs','Note','Nudged','Orientation','Best_Hit_Bitscore','Best_Hit_ARO',
                        'Predicted_DNA','Predicted_Protein','CARD_Protein_Sequence',
                        'ID','Percentage Length of Reference Sequence','Model_ID'])


    def parse_depcov_file(self):
        self.data = pd.read_csv(self.depcov,sep='\t',encoding='utf-8',header=None)
        self.data.columns = ["Contig","Depth","Coverage"]
        self.data["Depth"] = self.data["Depth"].round(2)
        self.data["Coverage"] = self.data["Coverage"].round(2)


    def parse_ref_file(self):
        ARO_list = set(self.df['ARO'].unique())
        self.ref = pd.read_csv(self.ref,sep='\t',encoding='utf-8',index_col=False)
        filter_ref = self.ref[self.ref['ARO'].isin(ARO_list)]
        print(self.ref)

        logging.info(f"[Output]: {self.outreffa}")
        with open(self.outreffa,'w',encoding='utf-8') as w:
            for index, row in filter_ref.iterrows():
                w.write(f">ARO:{row['ARO']}|Name:{row['Name']}|NCBI:{row['NCBI']}\n")
                w.write(f"{row['seq']}\n")
        

    def run(self):
        """输出"""
        self.parse_drug_anno()
        self.parse_depcov_file()
        self.parse_ref_file()

        out_df = pd.merge(self.df,self.data,on='Contig',how='inner')
        out_df.to_csv(self.table,sep='\t',encoding='utf-8',index=False)
        logging.info(f"[Output]: {self.table}")        



@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-1','--file', required=True,type=click.Path(exists=True),help="输入rgi.txt文件")
@click.option('-2','--depcov', required=True,type=click.Path(exists=True),help="输入深度覆盖度文件")
@click.option('--ref', required=True,type=click.Path(exists=True),help="参考序列tab")
def main(file,depcov,ref):
    project = Rgi(
        file = file,
        depcov = depcov,
        ref = ref
    )
    project.run()


if __name__ == '__main__':
    main()  