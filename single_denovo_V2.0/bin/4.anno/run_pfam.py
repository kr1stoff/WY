#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/22 16:40
import os
import logging
import click
import sys

logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


class Pfam:
    def __init__(self,**kwargs):
        self.out = kwargs['out']
        self.fa = kwargs['fasta']
        self.soft = kwargs['soft']
        self.db = kwargs['database']


    def make_dir(self):
        """makedir result dir"""
        os.makedirs(self.out,exist_ok=True)


    def check_file(self):
        """Check input path"""
        check_list = [self.fa, self.soft, self.db]
        for file in check_list:
            if not os.path.exists(os.path.abspath(file)):
                logging.error(f"No such file/dir [ {file} ]")
                exit(-1)


    def pfam_sh(self):
        """Write run pfam shell"""
        sh = os.path.join(self.out,'work.sh')
        with open(sh ,'w' ,encoding='utf-8') as w:
            w.write(f"{self.soft} -fasta {self.fa} -dir {self.db} -outfile {self.out}/pfam.out \n")
            w.write(f"grep -v '^#' {self.out}/pfam.out |sed '1d' >{self.out}/pfam.txt \n")
            w.write(f"sed -i '1i seq id\\talignment start\\talignment end\\tenvelope start\\tenvelope end\\thmm acc\\thmm name\\ttype\\thmm start\\thmm end\\thmm length\\tbit score\\tE-value\\tsignificance\\tclan' {self.out}/pfam.txt \n") 
        os.system(f"bash {sh}")
        logging.info(f"Output file [ {self.out}/pfam.txt ]")


#参数
p_soft = "/home/earthtest/miniconda3/envs/denovo/bin/pfam_scan.pl"
p_db = "/sdbb/share/database/Pfam"
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--fasta',required=True,type=click.Path(),help="fasta file")
@click.option('-o', '--out',required=True,type=click.Path(),help="Output dir")
@click.option('-p', '--soft',type=click.Path(),show_default=True,default=p_soft,help="Path of pfam_scan.pl ")
@click.option('-db', '--database',type=click.Path(),show_default=True,default=p_db,help="Path of Pfam database ")

def main(fasta,out,soft,database):
    """
    This pipeline need in denovo (conda'env)\n
    **Tips: Pipeline 's run Time >30min
    """

    project = Pfam(
                    fasta = fasta,
                    out = out,
                    soft = soft,
                    database = database
                    )
    project.make_dir()
    project.check_file()
    project.pfam_sh()



if __name__ == "__main__":
    main()




