#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/23 15:00
import os
import yaml
import sys
import click


# read YAML
PATH = os.path.dirname(os.path.abspath(__file__))  
DIR=os.path.dirname(PATH)
BIN = os.path.join(DIR, "bin") 
sys.path.append(DIR)
from lib.pipe import Pipe
import lib.common as common
cpu,th = common.get_cfg()
params = common.config()   # 读取配置文件


class VIRUS(Pipe):
    def __init__(self,run_flag=False,**kwargs):
        pipeline_dic = yaml.safe_load(open(kwargs['pipeline_yaml'] ))
        self.__dict__.update(pipeline_dic)
        self.run_flag = run_flag

    def make_dir(self):
        os.makedirs(self.p_sh ,exist_ok=True)
        os.makedirs(self.p_ass ,exist_ok=True)


    def align(self):
        sh  = os.path.join(self.p_sh,"2.ass.sh")
        cmd=[]
        cmd.append(f"""#Virus consensus
cd {self.p_ass}
source {params.ACTIVATE} denovo
{params.bwa} index -a bwtsw {self.ref}
""")
        if self.PE_flag:
            cmd.append(f"{params.bwa} mem -t {th} -Y {self.ref} {self.fq1} {self.fq2} \
                       |{params.samtools} view -@ 24 -hS -bF 12 \
                       |{params.samtools} sort -@ 24 -o {self.p_ass}/{self.id}_bf12_sort.bam")
        else:
            cmd.append(f"{params.bwa} mem -t {th} -Y {self.ref} {self.fq1} \
                       |{params.samtools} view -@ 24 -hS -bF 12 \
                       |{params.samtools} sort -@ 24 -o {self.p_ass}/{self.id}_bf12_sort.bam")

        cmd.append(f"""
{params.samtools} depth -a {self.p_ass}/{self.id}_bf12_sort.bam >{self.p_ass}/bf12_sort.depth
{params.samtools} mpileup -aa -A -d 0 -Q 0 {self.p_ass}/{self.id}_bf12_sort.bam | ivar consensus -p {self.id}_consensus
{params.samtools} mpileup -aa -A -d 0 -B -Q 0 --reference {self.ref} {self.p_ass}/{self.id}_bf12_sort.bam |ivar variants -p variants -r {self.ref}

#Depth dic
python3 {BIN}/2.ass_qc/depth_pic.py {self.p_ass}/bf12_sort.depth {self.p_ass} {self.ref}
python3 {BIN}/2.ass_qc/ivar_var.py variants.tsv {self.p_ass}
""")
        common.cmd2shell(cmd,sh)
        return sh


    def run(self):
        self.make_dir()
        sh = self.align()
        if not self.run_flag:
            common.prun((sh ,self.log))



# == Option ================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--pipeline_yaml',type=click.Path(),help="分析目录下的pipeline.yaml文件")
@click.option('--norun',is_flag=True,help="是否只生成命令，不运行程序,默认运行")
def main(pipeline_yaml,norun):
    project = VIRUS(
        pipeline_yaml = pipeline_yaml,
        run_flag = norun
    )
    project.run()



if __name__ == "__main__":
    main()