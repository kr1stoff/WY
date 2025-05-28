#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/22 16:40
import os
import logging
import re
import click
import sys

logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


class PhiSpy:
    def __init__(self,res):
        self.res = res
        os.makedirs(self.res,exist_ok=True)

    def parse_gbk(self,gbkFile):
        """
        PhiSpy don't support long name in gbk file.
        """
        with open(gbkFile ,'r',encoding='utf-8') as f ,open(f"{self.res}/PhiSpy.gbk" ,'w' ,encoding='utf-8') as w:
            for line in f:
                if not line: continue
                if re.search('^LOCUS' ,line):
                    li = re.split('\s+',line)
                    id = li[1]
                    tail = li[-3] +' '+ li[-2]
                    tmp = re.findall(r"(.*)_length_(\d+)_" ,id)[0]
                    short_id = tmp[0]
                    length = tmp[1]
                    w.write(f"LOCUS       {short_id} {length} bp   {tail}  ")
                else:
                    w.write(line)


    def phispy_sh(self,soft):
        """Write run phispy shell"""
        sh = os.path.join(self.res,'work.sh')
        with open(sh ,'w' ,encoding='utf-8') as w:
            w.write(f"{soft} {self.res}/PhiSpy.gbk -u 1000 --output_choice 511 -o {self.res} \n")  #/home/earthtest/miniconda3/envs/denovo/bin/PhiSpy.py
            w.write( f"sort -Vk 1 {self.res}/prophage_coordinates.tsv |cut -f 1-8,11|sed -e '1i 前噬菌体ID\\t组装序列ID\\tprophage start location\\tprophage stop location\\tstart of attL\\tend of attL\\tstart of attR\\tend of attR\\t说明' > {self.res}/prophage_coordinates.txt \n")
        os.system(f"bash {sh}")
        logging.info(f"Output file [{self.res}/prophage_coordinates.txt]")


#参数
p_soft = "/home/earthtest/miniconda3/envs/denovo/bin/PhiSpy.py"
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gbk',required=True,type=click.Path(),help="Prokka gbk file")
@click.option('-o', '--out',required=True,type=click.Path(),help="Output dir")
@click.option('-p', '--path',type=click.Path(),show_default=True,default=p_soft,help="PhiSpy software path")

def main(gbk,out,path):
    """
    1.This pipeline need in denovo (conda'env)\n
    2.This pipeline need contig_headers such as [NODE_1_length_3451_cov_22]
    """

    project = PhiSpy(out)
    project.parse_gbk(gbk)
    project.phispy_sh(p_soft)



if __name__ == "__main__":
    main()




