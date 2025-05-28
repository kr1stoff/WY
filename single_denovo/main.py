#!/home/earthtest/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2023/5/10
import click
import os
import sys
import yaml
PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(PATH))

import lib.common as common
from lib.pipe import Pipe
from modules.SeqQc import SeqQc
from modules.Assemble import Assemble
from modules.AssQC import AssembleQC
from modules.Predict import Predict
from modules.Anno import Anno
from modules.Virus import VIRUS
from modules.Report import Report

params = common.config()   # 读取配置文件

class Denovo(Pipe):
    def __init__(self,**kwargs):
        id = kwargs['id']
        ref = kwargs['ref']
        solar_yaml = os.path.abspath(kwargs['solar_yaml'])
        
        # 读取 lib/Pipe 模块，生成该样本的 pipeline.yaml 
        project = Pipe(solar_yaml = solar_yaml, id = id, ref = ref)
        project.run()

        # 读取该样本的 pipeline.yaml ，转化为 self 变量
        self.pipeline_yaml = os.path.abspath(kwargs['pipeline_yaml'])
        pipeline_dic = yaml.safe_load(open(self.pipeline_yaml))
        self.__dict__.update(pipeline_dic)


    def run(self):
        # 细菌/真菌
        if self.kind == 'bacteria' or self.kind == 'fungi':
            project1 = SeqQc(pipeline_yaml = self.pipeline_yaml)
            project2 = Assemble(pipeline_yaml = self.pipeline_yaml)     # 细菌进行截取 9M 数据，真菌/病毒不截
            project3 = AssembleQC(pipeline_yaml = self.pipeline_yaml)
            project4 = Predict(pipeline_yaml = self.pipeline_yaml)      # 细菌 真菌，预测软件分开
            project5 = Anno(pipeline_yaml = self.pipeline_yaml)
            project6 = Report(pipeline_yaml = self.pipeline_yaml)
            project1.run()
            project2.run()
            project3.run()    
            project4.run() 
            project5.run() 
            project6.run()

        # 病毒有参拼接
        elif self.kind == 'virus'and not self.if_denovo:
            project1 = SeqQc(pipeline_yaml = self.pipeline_yaml)
            project2 = VIRUS(pipeline_yaml = self.pipeline_yaml)
            project3 = Report(pipeline_yaml = self.pipeline_yaml)
            project1.run()
            project2.run()
            project3.run()

        # 病毒无参组装
        elif self.kind == 'virus' and self.if_denovo:
            project1 = SeqQc(pipeline_yaml = self.pipeline_yaml)
            project2 = Assemble(pipeline_yaml = self.pipeline_yaml)
            project3 = AssembleQC(pipeline_yaml = self.pipeline_yaml)
            project4 = Report(pipeline_yaml = self.pipeline_yaml)
            project1.run()
            project2.run()
            project3.run()
            project4.run()



####参数#########################################################################################
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i','--id', required=True,type=click.STRING,help="样本名称")
@click.option('-p','--pipeline_yaml',type=click.Path(),help="输入分析目录下的pipeline.yaml文件")
@click.option('-s','--solar_yaml',type=click.Path(),help="solar yaml")
@click.option('-r','--ref',default='None',type=click.STRING,help="参考序列，默认为None")

def main(id,pipeline_yaml,solar_yaml,ref):
    """创建单菌流程中，需要的分析目录"""
    project = Denovo(
        id = id,
        pipeline_yaml = pipeline_yaml,
        solar_yaml = solar_yaml,
        ref = ref
        )
    project.run()



if __name__ == "__main__":
    main()