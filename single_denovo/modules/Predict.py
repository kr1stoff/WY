#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/23 15:00
import os
import sys
import click
import yaml

# read YAML
PATH = os.path.dirname(os.path.abspath(__file__))  
DIR=os.path.dirname(PATH)
BIN = os.path.join(DIR, "bin") 
sys.path.append(DIR)
from lib.pipe import Pipe
import lib.common as common
cpu,th = common.get_cfg()
params = common.config()   # 读取配置文件

class Predict(Pipe):
    def __init__(self,**kwargs):
        pipeline_dic = yaml.safe_load(open(kwargs['pipeline_yaml'] ))
        self.__dict__.update(pipeline_dic)
        self.infa = self.filter_contig   #输入为过滤后的fa

    def creat_env(self):
        os.makedirs(self.p_predict,exist_ok=True)
        os.makedirs(self.p_island,exist_ok=True)
        # os.makedirs(self.p_phispy,exist_ok=True)


    def bacteria_predict(self):
        cmd = []
        sh = os.path.join(self.p_sh,"4.precict.sh")
        cmd.append(f"""#prokka
source {params.ACTIVATE} denovo
cd {self.p_predict}
time prokka {self.infa} --prefix predict --cpus {int(cpu*3/4)} --outdir {self.p_predict} --kingdom Bacteria --force --addgenes --quiet

if [ -f "{self.p_predict}/predict.gff" ]; then
    {params.python3} {BIN}/3.pre/chinese_prokka.py -i {self.p_predict}/predict.txt -o {self.p_predict}/predict_kind.txt

    # bed
    grep -v '#' {self.p_predict}/predict.gff|grep 'CDS' |grep 'ID='|cut -f 1,4,5,9|awk -F';' '{{print $1}}'|sed 's/ID=//g' >{self.p_predict}/predict.bed
    bedtools coverage -mean -a {self.p_predict}/predict.bed -b {self.sort_bam} >{self.p_predict}/bed_bam.merge.depth
    bedtools coverage -a {self.p_predict}/predict.bed -b {self.sort_bam} >{self.p_predict}/bed_bam.merge.cov
    paste {self.p_predict}/bed_bam.merge.depth {self.p_predict}/bed_bam.merge.cov |cut -f 4,5,13 >{self.merge_file}

    # gene_length
    {params.b_denovo}/python3 {BIN}/2.ass_qc/ass_fa_stat.py -fa {self.p_predict}/predict.ffn -p predict
fi

    # island
    time {params.python3} {BIN}/3.pre/multi_Island.py -i {self.p_predict}/predict.gbk -o {self.p_island} >{self.p_island}/island.log 2>&1
    #phispy
    time {params.python3} {BIN}/3.pre/run_phispy.py -i {self.p_predict}/predict.gbk -o {self.p_phispy}
""")

        common.cmd2shell(cmd,sh)
        return sh

    
    def fungi_predict(self):
        cmd = []
        sh = os.path.join(self.p_sh,"4.precict.sh")
        cmd.append(f"""#genmark
source {params.ACTIVATE} denovo
cd {self.p_predict}
time {params.genemark} --ES --fungus --cores {int(cpu*3/4)} --sequence {self.infa}  >{self.p_predict}/pre.log

if [ -f "{self.p_predict}/genemark.gtf" ]; then
    echo "gtf ==> gff3 ==> cds.fa ==> cds.faa"
    gffread {self.p_predict}/genemark.gtf -o ->{self.p_predict}/predict.gff3
    {params.perl} {BIN}/tools/getGene.pl {self.p_predict}/predict.gff3 {self.filter_contig} -type mrna >{self.p_predict}/predict.fa
    {params.perl} {BIN}/tools/cds2aa.pl {self.p_predict}/predict.fa >{self.p_predict}/predict.faa

    grep 'CDS' {self.p_predict}/predict.gff3|grep 'Parent='|cut -f 1,4,5,9|sed -e 's/Parent=//g' >{self.p_predict}/predict.gff3.bed
    bedtools coverage -mean -a {self.p_predict}/predict.gff3.bed -b {self.bam} >{self.p_predict}/bed_bam.merge.depth
    bedtools coverage -a {self.p_predict}/predict.gff3.bed -b {self.bam} >{self.p_predict}/bed_bam.merge.cov
    paste {self.p_predict}/bed_bam.merge.depth {self.p_predict}/bed_bam.merge.cov |cut -f 4,5,13 >{self.merge_file}

    {params.python3} {BIN}/3.pre/gene_mark_count.py {self.p_predict}/predict.gtf {self.p_predict}/predict_kind.txt
    {params.b_denovo}/python3 {BIN}/2.ass_qc/ass_fa_stat.py -fa {self.p_predict}/predict.fa -p predict
fi
""")

        common.cmd2shell(cmd,sh)
        return sh
    
    def run(self):
        self.creat_env()
        if self.kind == 'bacteria': sh = self.bacteria_predict()
        if self.kind == 'fungi':    sh = self.fungi_predict()
        common.prun((sh, self.log))


# ==== Option =============================================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--pipeline_yaml',type=click.Path(),help="分析目录下的pipeline.yaml文件")

def main(pipeline_yaml):
    project = Predict(
            pipeline_yaml = pipeline_yaml,
            )
    project.run()




if __name__ == "__main__":
    main()
