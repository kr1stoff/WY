#!/usr/bin/env python
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

class AssembleQC(Pipe):
    def __init__(self,**kwargs):
        pipeline_dic = yaml.safe_load(open(kwargs['pipeline_yaml'] ))
        self.__dict__.update(pipeline_dic)
        self.infa = self.filter_contig   #输入为过滤后的fa


    def creat_env(self):
        os.makedirs(self.p_quast ,exist_ok=True)
        os.makedirs(self.p_checkm ,exist_ok=True)
        os.makedirs(self.p_depth ,exist_ok=True)
        if self.ref:
            os.makedirs(self.p_ref ,exist_ok=True)


    def quast_checkm(self):
        sh = os.path.join(self.p_sh,"3.1quast_checkm.sh")
        cmd = []
        cmd.append(f"""
source {params.ACTIVATE} denovo
python3 {params.quast} {self.infa} --plots-format png -o {self.p_quast} --silent -t {int(th/2)}
checkm lineage_wf -t {int(th*3/4)} -q --tab_table -f {self.p_checkm}/checkm.txt -x fasta {self.p_spades} {self.p_checkm}  >{self.p_checkm}/checkm.o 2>&1 
{params.python3} {BIN}/2.ass_qc/parse_checkm_result.py -i {self.p_checkm}/storage/bin_stats_ext.tsv -o {self.p_checkm}/check.txt 
""")
        common.cmd2shell(cmd,sh)
        return sh



    def depth(self):
        sh = os.path.join(self.p_sh,"3.2depth_stat.sh")
        cmd = []
        
        cmd.append(f"{params.bwa} index -a bwtsw {self.infa}")
        cmd.append(f"{params.samtools} faidx {self.infa}")

        if self.PE_flag:  #PE
            cmd.append(f"mode=PE")
            cmd.append(f"{params.bwa} mem -t {int(th*3/4)} {self.infa} {self.fq1} {self.fq2} \
                        |{params.samtools} view -@ 24 -bS \
                        |{params.samtools} sort -@ 24 >{self.sort_bam}")
        else:   #SE
            cmd.append(f"mode=SE")
            cmd.append(f"{params.bwa} mem -t {int(th*3/4)} {self.infa} {self.fq1} \
                        |{params.samtools} view -@ 24 -bS \
                        |{params.samtools} sort -@ 24 >{self.sort_bam}")

        cmd.append(f"""#samtools
{params.samtools} depth -a {self.sort_bam} >{self.p_depth}/sort.bam.depth

# Insert size
if [ "$mode" == 'PE' ];then
    {params.picard} CollectInsertSizeMetrics --Histogram_FILE {self.p_depth}/picard_insertsize.pdf -I {self.sort_bam} -O {self.p_depth}/insertsize.txt
    {params.Rscript} {BIN}/2.ass_qc/insertsize.R {self.p_depth}/insertsize.txt {self.p_depth}
fi

# GC-depth png
{params.python3} {BIN}/2.ass_qc/depth_base_stat.py -l 2000 -g {self.infa} -d {self.p_depth}/sort.bam.depth -s {self.p_depth}/depth_base.stat
{params.Rscript} {BIN}/2.ass_qc/depth_GC_plot.r -i {self.p_depth}/depth_base.stat -o {self.p_depth}/depth_base.stat.depth_GC
{params.convert} {self.p_depth}/depth_base.stat.depth_GC.pdf {self.p_depth}/depth_base.stat.depth_GC.png

# Uniformity
{params.python3} {BIN}/2.ass_qc/count_depth.py -i {self.p_depth}/sort.bam.depth -ref {self.infa} -o {self.p_depth}/uniformity.txt -a {self.p_depth}/uniformity_all.txt
""")
        common.cmd2shell(cmd,sh)
        return sh



    def align_ref(self):
        sh = os.path.join(self.p_sh,"3.3align_ref.sh")
        cmd = []
        if self.ref:
            ref_name = os.path.basename(self.ref)
            cmd.append(f"""
cd {self.p_ref}
cp {self.ref} {self.p_ref}
{params.makeblastdb} -in {self.p_ref}/{ref_name} -out ref -dbtype nucl
{params.blastn} -num_threads {int(th/2)} -evalue 1e-10 -query {self.infa} -db {self.p_ref}/ref -outfmt '{params.blast_opt}' -out {self.p_ref}/ref_blastn.txt
{params.python3} {BIN}/tools/cal_blast.py -i {self.p_ref}/ref_blastn.txt --mismatch_rate 0.1 --identity_rate 0.8
""")
        common.cmd2shell(cmd,sh)
        return sh



    def run(self):
        self.creat_env()
        sh1 = self.quast_checkm()
        sh2 = self.depth()
        sh3 = self.align_ref()
        shlist=[sh1,sh2,sh3]
        common.mul_pool(shlist, self.log)


# Option
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--pipeline_yaml',type=click.Path(),help="分析目录下的pipeline.yaml文件")
def main(pipeline_yaml):
    project = AssembleQC(
        pipeline_yaml = pipeline_yaml
            )
    project.run()



if __name__ == "__main__":
    main()