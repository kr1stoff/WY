#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/23 15:00
import os
import yaml
import sys
import click

PATH = os.path.dirname(os.path.abspath(__file__))  
DIR=os.path.dirname(PATH)
BIN = os.path.join(DIR, "bin") 
sys.path.append(DIR)
import lib.common as common
cpu,th =  common.get_cfg()
params = common.config()   # 读取配置文件
from lib.pipe import Pipe


class SeqQc(Pipe):
    def __init__(self,**kwargs):
        pipeline_dic = yaml.safe_load(open(kwargs['pipeline_yaml'] ))
        self.__dict__.update(pipeline_dic)

    def creat_env(self):
        os.makedirs(self.p_fastqc,exist_ok=True)
        os.makedirs(self.p_fastp,exist_ok=True)
        os.makedirs(self.p_filter_fastqc,exist_ok=True)

    def fastqc(self):
        """ run fastp """
        sh = os.path.join(self.p_sh,"1.1fastqc.sh")
        cmd = []
        
        if self.PE_flag: #PE
            l_fq1 ,l_fq2 = common.pe_fq_link(self.id,self.fq1,self.fq2,self.p_fastqc) 
            cmd.append(f"{params.fastqc} -q -t {th} {l_fq1} {l_fq2} -o {self.p_fastqc}")
            cmd.append(f"unzip -qn {self.p_fastqc}/*1_fastqc.zip -d {self.p_fastqc}")
            cmd.append(f"unzip -qn {self.p_fastqc}/*2_fastqc.zip -d {self.p_fastqc}")
        else: #SE
            l_fq = common.se_fq_link(self.id,self.fq1,self.p_fastqc)
            cmd.append(f"{params.fastqc} -q -t {th} {l_fq}  -o {self.p_fastqc}")
            cmd.append(f"unzip -qn {self.p_fastqc}/*_fastqc.zip -d {self.p_fastqc}")

        common.cmd2shell(cmd,sh)
        return sh


    def fastp(self):
        """ run fastp """
        sh = os.path.join(self.p_sh,"1.2fastp.sh")
        cmd=[]
        if self.PE_flag: #PE
            cmd.append(f"{params.fastp} -n 1 -l 50 -y -w 16 -i {self.fq1} -I {self.fq2} -o {self.clean_fq1} -O {self.clean_fq2} -j {self.p_fastp}/fastp.json  -h {self.p_fastp}/fastp.html ")
            cmd.append(f"{params.python3} {BIN}/1.qc/parse_fastp.json.py -id {self.id} -i {self.p_fastp}/fastp.json -m PE ")

        else: #SE
            cmd.append(f"{params.fastp} -n 1 -l 50 -y -w 16 -i {self.fq1} -j {self.p_fastp}/fastp.json -h {self.p_fastp}/fastp.html -o {self.clean_fq1} ")
            cmd.append(f"{params.python3} {BIN}/1.qc/parse_fastp.json.py -id {self.id} -i {self.p_fastp}/fastp.json -m SE ")

        common.cmd2shell(cmd,sh)
        return sh

    def filter_fastqc(self):
        """ run fastqc """
        sh = os.path.join(self.p_sh,"1.3filter_fastqc.sh")
        cmd = []
        
        if self.PE_flag: #PE
            cmd.append(f"{params.fastqc} -q -t {th} {self.clean_fq1} {self.clean_fq2} -o {self.p_filter_fastqc}")
            cmd.append(f"unzip -qn {self.p_filter_fastqc}/*_1.clean_fastqc.zip -d {self.p_filter_fastqc}")
            cmd.append(f"unzip -qn {self.p_filter_fastqc}/*_2.clean_fastqc.zip -d {self.p_filter_fastqc}")
        else: #SE
            cmd.append(f"{params.fastqc} -q -t {th} {self.clean_fq1} -o {self.p_filter_fastqc}")
            cmd.append(f"unzip -qn {self.p_filter_fastqc}/*_fastqc.zip -d {self.p_filter_fastqc}")

        common.cmd2shell(cmd,sh)
        return sh


    def run(self):
        """运行"""
        self.creat_env()
        sh1 = self.fastqc()
        sh2 = self.fastp()
        shlist1=[sh1,sh2]
        common.mul_pool(shlist1, self.log)

        sh4 = self.filter_fastqc()
        common.prun((sh4, self.log))



# # === Option ===============================================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--pipeline_yaml',type=click.Path(),help="单跑，仅输入分析目录下的pipeline.yaml文件")

def main(pipeline_yaml):
    project = SeqQc(
        pipeline_yaml = pipeline_yaml
    )
    project.run()


if __name__ == "__main__":
    main()
