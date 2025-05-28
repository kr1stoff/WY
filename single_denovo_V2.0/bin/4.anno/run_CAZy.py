#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/22 16:40
import os
import logging
import click
import sys

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
p_des = "/sdbb/share/database/CAZY/CAZyDB.07302020.fam-activities-hmm.txt"

class CAZy:
    def __init__(self,**kwargs):
        self.out = kwargs['out']
        self.faa = kwargs['faa']
        self.soft = kwargs['soft']
        self.db = kwargs['database']
        #self.cpu = kwargs['cpu']
        self.tool = kwargs['tool']


    def make_dir(self):
        """makedir result dir"""
        os.makedirs(self.out,exist_ok=True)


    def check_file(self):
        """Check input path"""
        check_list = [self.faa, self.soft, self.db, p_des]
        for file in check_list:
            if not os.path.exists(os.path.abspath(file)):
                logging.error(f"No such file/dir [ {file} ]")
                exit(-1)


    def CAZy_sh(self):
        """Write run CAZy shell"""
        sh = os.path.join(self.out,'work.sh')
        with open(sh ,'w' ,encoding='utf-8') as w:
            w.write(f"{self.soft} {self.faa} protein --db_dir {self.db}  --out_dir {self.out} --tools {self.tool} \n")
            w.write(f"csvtk join -t -L -f 1  {self.out}/hmmer.out {p_des} |csvtk round -t -f 'Coverage' -n 2 |csvtk mutate -t -f 1 -p '^(\D+)' -n 'kind' >{self.out}/CAZY.txt \n")

        os.system(f"bash {sh}")



#参数
p_soft = "/home/earthtest/miniconda3/envs/denovo/bin/run_dbcan"
p_db = "/sdbb/share/database/CAZY"
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--faa',required=True,type=click.Path(),help="protein faa")
@click.option('-o', '--out',required=True,type=click.Path(),help="Output dir")
#@click.option('--cpu',type=click.INT,show_default=True,default=16,help="Numbers of use cpu cores ")
@click.option('-p', '--soft',type=click.Path(),show_default=True,default=p_soft,help="Path of run_dbcan ")
@click.option('-db', '--database',type=click.Path(),show_default=True,default=p_db,help="Path of CAZy database ")
@click.option('--tool',type=click.STRING,show_default=True,default='hmmer',help="choose predict tools {hmmer,diamond,eCAMI,all ")

def main(faa,out,soft,database,tool):
    """
    This pipeline need in denovo (conda'env)\n
    **Tips: Pipeline 's run Time ~15min \n
    No clear change even if change cpu numbers
    """

    project = CAZy( faa = faa,
                    out = out,
                    soft = soft,
                    database = database,
                    tool = tool
                    )
    project.make_dir()
    project.check_file()
    project.CAZy_sh()



if __name__ == "__main__":
    main()