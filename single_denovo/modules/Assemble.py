#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/23 15:00
import os
import yaml
import sys
import click
import sys
 
PATH = os.path.dirname(os.path.abspath(__file__))  
DIR=os.path.dirname(PATH)
BIN = os.path.join(DIR, "bin") 
sys.path.append(DIR)
from lib.pipe import Pipe
import lib.common as common
cpu,th = common.get_cfg()
params = common.config()   # 读取配置文件

class Assemble(Pipe):
    def __init__(self,**kwargs):
        pipeline_dic = yaml.safe_load(open(kwargs['pipeline_yaml'] ))
        self.__dict__.update(pipeline_dic)


    def creat_env(self):
        os.makedirs(self.p_deal,exist_ok=True)
        os.makedirs(self.p_spades,exist_ok=True)


    def Advanced_Option(self):
        """ 去污染/截断reads """
        sh = os.path.join(self.p_sh,"2.1Advanced_Option.sh")
        cmd = []
        self.in_fq1 ,self.in_fq2 = self.clean_fq1 ,self.clean_fq2   # 定义组装输入fq

        if self.PE_flag:  #PE
            cmd.append(f"{params.python3} {BIN}/1.qc/depollute/depollute_by_kraken2.py -p {self.p_deal}/depollute --config {BIN}/1.qc/depollute/depollute_by_kraken2.yaml {self.in_fq1} {self.in_fq2}")

        else: #SE
            cmd.append(f"{params.python3} {BIN}/1.qc/depollute/depollute_by_kraken2.py -p {self.p_deal}/depollute --config {BIN}/1.qc/depollute/depollute_by_kraken2.yaml {self.clean_fq1}")
        cmd.append(f"{params.python3} {BIN}/1.qc/parse_kraken.report.py -i {self.p_deal}/temp_depollute/depollute.report ")
        
        # 去污染
        if self.if_pollute:
            self.in_fq1 ,self.in_fq2 = self.depollute_fq1 ,self.depollute_fq2   # 定义组装输入fq


        # 细菌进行截reads
        if self.kind == 'bacteria':    
            if self.PE_flag:  #PE
                cmd.append(f"{params.seqtk} sample -s 11 {self.in_fq1} {params.cut_reads} >{self.cutfq1} ")
                cmd.append(f"{params.seqtk} sample -s 11 {self.in_fq2} {params.cut_reads} >{self.cutfq2} ")
                cmd.append(f"{params.iTools} Fqtools stat -InFq {self.cutfq1} -InFq {self.cutfq2} -OutStat {self.p_deal}/cutfq.stat")
                cmd.append(f"{params.python3} {BIN}/1.qc/parse.itools.stat.py -i {self.p_deal}/cutfq.stat -m PE")

            else:   #SE
                cmd.append(f"{params.seqtk} sample -s 11 {self.in_fq1} {params.cut_reads} >{self.cutfq1}")
                cmd.append(f"{params.iTools} Fqtools stat -InFq {self.cutfq1}  -OutStat {self.p_deal}/cutfq.stat")
                cmd.append(f"{params.python3} {BIN}/1.qc/parse.itools.stat.py -i {self.p_deal}/cutfq.stat -m SE")
            self.in_fq1 ,self.in_fq2 = self.cutfq1 ,self.cutfq2     # 定义组装输入fq

        common.cmd2shell(cmd,sh)
        return sh



    def spades(self):
        sh = os.path.join(self.p_sh,"2.2spades.sh")
        cmd=[]
        cmd.append(f"source {params.ACTIVATE} denovo")
        if self.PE_flag:  #PE
            cmd.append(f"time spades.py -t {int(th*3/4)} --careful -o {self.p_spades} -1 {self.in_fq1} -2 {self.in_fq2} >{self.p_spades}/spades.o ")
        else: #SE
            cmd.append(f"time spades.py -t {int(th*3/4)} --careful -o {self.p_spades} -s {self.in_fq1} >{self.p_spades}/spades.o ")
        cmd.append(f"""
{params.python3} {BIN}/2.ass_qc/ass_fa.filter.py -fa {self.p_spades}/scaffolds.fasta -p {self.id}_scaffolds -l {params.minlen}
{params.python3} {BIN}/2.ass_qc/ass_fa_stat.py -fa {self.p_spades}/{self.id}_scaffolds.fasta -p {self.id}_scaffolds
""")
        common.cmd2shell(cmd,sh)
        return sh



    def run(self):
        self.creat_env()
        sh1 = self.Advanced_Option()
        sh2 = self.spades()
        common.prun((sh1, self.log))
        common.prun((sh2, self.log))




# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--pipeline_yaml',type=click.Path(),help="分析目录下的pipeline.yaml文件")
def main(pipeline_yaml):
    project = Assemble(
        pipeline_yaml = pipeline_yaml
    )
    project.run()



if __name__ == "__main__":
    main()