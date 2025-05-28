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
params = common.config()   # 读取配置文件

class Report(Pipe):
    def __init__(self,run_flag=False,**kwargs):
        pipeline_dic = yaml.safe_load(open(kwargs['pipeline_yaml'] ))
        self.__dict__.update(pipeline_dic)
        self.run_flag = run_flag

    def make_dir(self):
        if self.kind != 'virus':
            for var in ['1.qc','2.ass','3.pre','4.anno']:   #report dir
                tmp = os.path.join(self.p_report,var)
                os.makedirs(tmp,exist_ok=True)
        else:
            for var in ['1.qc','2.ass']:   #report dir
                tmp = os.path.join(self.p_report,var)
                os.makedirs(tmp,exist_ok=True)


    def run_report(self):
        sh  = os.path.join(self.p_sh,"6.report.sh")
        cmd = []
        cmd.append(f"""
cd {self.out}

# 1.qc
cp -r {self.p_fastp}/{self.id}.basic.stat.txt {self.p_report}/1.qc

if [ -f "{self.p_fastqc}/{self.id}_2_fastqc.html" ]; then
    cp -r {self.p_fastqc}/{self.id}_1_fastqc {self.p_report}/1.qc
    cp -r {self.p_fastqc}/{self.id}_1_fastqc.html {self.p_report}/1.qc
    cp -r {self.p_filter_fastqc}/{self.id}_1.clean_fastqc {self.p_report}/1.qc
    cp -r {self.p_filter_fastqc}/{self.id}_1.clean_fastqc.html {self.p_report}/1.qc


    cp -r {self.p_fastqc}/{self.id}_2_fastqc {self.p_report}/1.qc
    cp -r {self.p_fastqc}/{self.id}_2_fastqc.html {self.p_report}/1.qc
    cp -r {self.p_filter_fastqc}/{self.id}_2.clean_fastqc {self.p_report}/1.qc
    cp -r {self.p_filter_fastqc}/{self.id}_2.clean_fastqc.html {self.p_report}/1.qc

else
    cp -r {self.p_fastqc}/{self.id}_fastqc {self.p_report}/1.qc
    cp -r {self.p_fastqc}/{self.id}_fastqc.html {self.p_report}/1.qc
    cp -r {self.p_filter_fastqc}/{self.id}.clean_fastqc {self.p_report}/1.qc
    cp -r {self.p_filter_fastqc}/{self.id}.clean_fastqc.html {self.p_report}/1.qc

fi
""")
 
        # ========================================================================
        if self.kind == 'bacteria' or self.kind == 'fungi':
            cmd.append(f"""
# 2.ass
cp -r {self.p_deal}/cutfq.stat.txt  {self.p_report}/1.qc
# depollute
if  [ -f "{self.p_deal}/temp_depollute/depollute.report" ];then 
    cp -r {self.p_deal}/temp_depollute/abundance.txt  {self.p_report}/1.qc
fi

cp -r {self.p_spades}/{self.id}_scaffolds_fa.stat.txt {self.p_report}/2.ass
cp -r {self.p_spades}/scaffolds.fasta    {self.p_report}/2.ass
cp -r {self.p_spades}/{self.id}_scaffolds.fasta    {self.p_report}/2.ass
cp -r {self.p_spades}/{self.id}_scaffolds.length.png    {self.p_report}/2.ass

cp -r {self.p_quast}            {self.p_report}/2.ass
cp -r {self.p_checkm}/check.txt {self.p_report}/2.ass
cp -r {self.p_depth}/depth_base.stat.depth_GC.png {self.p_report}/2.ass

if [ -f "{self.p_fastqc}/{self.id}_1_fastqc.html" ]; then
    cp -r {self.p_depth}/insertsize.png {self.p_report}/2.ass
fi

cp -r {self.p_depth}/depth_base.stat    {self.p_report}/2.ass
cp -r {self.p_depth}/uniformity.txt     {self.p_report}/2.ass
cp -r {self.p_depth}/uniformity_all.txt {self.p_report}/2.ass
cp -r {self.p_ref}/ref_cov_stat.txt {self.p_report}/2.ass
""")

            cmd.append(f"""
# 3.precict
cp -r {self.p_predict}/predict_kind.txt {self.p_report}/3.pre
cp -r {self.p_predict}/predict.length.png {self.p_report}/3.pre
cp -r {self.p_repeat}/INEs.txt {self.p_report}/3.pre
cp -r {self.p_repeat}/TandemRepeat.txt {self.p_report}/3.pre
""")

        # ========================================================================
        if self.kind == 'bacteria':
            cmd.append(f"""
cp -r {self.p_island}/island.txt {self.p_report}/3.pre
cp -r {self.p_phispy}/prophage_coordinates.txt {self.p_report}/3.pre
cp -r {self.p_phispy}/phage.fasta {self.p_report}/3.pre
""")

        # ========================================================================
        if self.kind == 'bacteria' or self.kind == 'fungi':
            cmd.append(f"""# 4.anno
cp -r {self.p_VFDB}/Virulence.txt  {self.p_report}/4.anno
cp -r {self.p_VFDB}/Virulence_gene.fa     {self.p_report}/4.anno
cp -r {self.p_VFDB}/Virulence_ref.fa {self.p_report}/4.anno

cp -r {self.p_drug}/CARD_anno.txt  {self.p_report}/4.anno
cp -r {self.p_drug}/CARD_contig.fa  {self.p_report}/4.anno
cp -r {self.p_drug}/CARD_ref.fa  {self.p_report}/4.anno
cp -r {self.p_drug}/Resfinder_anno.txt  {self.p_report}/4.anno
cp -r {self.p_drug}/ResFinder_Hit_in_genome_seq.fsa  {self.p_report}/4.anno

cp -r {self.p_eggmapper}/COG/all.COG.bar.png {self.p_report}/4.anno
cp -r {self.p_eggmapper}/COG/all.COG.class.xls {self.p_report}/4.anno
cp -r {self.p_eggmapper}/GO/GO_anno_stats.xls {self.p_report}/4.anno
cp -r {self.p_eggmapper}/GO/GO_anno.xls {self.p_report}/4.anno
cp -r {self.p_eggmapper}/GO/GO_anno_stats_level2.png {self.p_report}/4.anno
cp -r {self.p_eggmapper}/KEGG/KEGG_anno.txt {self.p_report}/4.anno
cp -r {self.p_eggmapper}/KEGG/KEGG_anno.stat.txt {self.p_report}/4.anno
cp -r {self.p_eggmapper}/KEGG/KEGG_anno_stats.png {self.p_report}/4.anno

cp -r {self.p_CAZy}/CAZY.txt {self.p_report}/4.anno
cp -r {self.p_CAZy}/CAZy_anno_stats.png {self.p_report}/4.anno
cp -r {self.p_swissprot}/swissprot_result.tsv {self.p_report}/4.anno

cp -r {self.p_plasmid}/plasmid.txt {self.p_report}/4.anno
""")

        # 病毒有参拼接 # ========================================================================
        elif self.kind == 'virus' and os.path.exists(f"{self.p_ass}/{self.id}_consensus.fa"):
            cmd.append(f"""
cp -r {self.p_ass}/{self.id}_bf12_sort.bam {self.p_report}/2.ass
cp -r {self.p_ass}/bed.depth.stat.png {self.p_report}/2.ass
cp -r {self.p_ass}/bam.stat.txt  {self.p_report}/2.ass
cp -r {self.p_ass}/variants.xlsx {self.p_report}/2.ass
cp -r {self.p_ass}/variants.txt  {self.p_report}/2.ass
cp -r {self.p_ass}/{self.id}_consensus.fa  {self.p_report}/2.ass
""")

        # 病毒无参拼接 # ========================================================================
        elif self.kind == 'virus' and os.path.exists(f"{self.p_spades}/{self.id}_scaffolds.fasta"):
            cmd.append(f"""
cp -r {self.p_spades}/{self.id}_scaffolds_fa.stat.txt {self.p_report}/2.ass
cp -r {self.p_spades}/scaffolds.fasta    {self.p_report}/2.ass
cp -r {self.p_spades}/{self.id}_scaffolds.fasta    {self.p_report}/2.ass
cp -r {self.p_spades}/{self.id}_scaffolds.length.png    {self.p_report}/2.ass
cp -r {self.p_quast}           {self.p_report}/2.ass
cp -r {self.p_depth}/depth_base.stat.depth_GC.png {self.p_report}/2.ass
if [ -f "{self.p_fastqc}/{self.id}_1_fastqc.html" ]; then
    cp -r {self.p_depth}/insertsize.png {self.p_report}/2.ass
fi

cp -r {self.p_depth}/depth_base.stat    {self.p_report}/2.ass
cp -r {self.p_depth}/uniformity.txt     {self.p_report}/2.ass
cp -r {self.p_depth}/uniformity_all.txt {self.p_report}/2.ass
cp -r {self.p_ref}/ref_cov_stat.txt     {self.p_report}/2.ass
""")

        # src # ========================================================================
        pipe_src = os.path.join(DIR,"lib","report","pipe_src")
        cmd.append(f"cp -r {params.p_GDHR}/src {self.p_report}")
        cmd.append(f"cp -r {pipe_src}/{self.kind}.png {self.p_report}/src/image")
        cmd.append(f"cp -r {pipe_src}/icon {self.p_report}/src/image")

        # 报告 # ========================================================================
        report_pl = os.path.join(DIR,"lib","report",f"{self.kind}_report.pl")
        cmd.append(f"{params.perl} {report_pl} {self.id} {self.out}")
        cmd.append(f"zip -qr {self.id}.zip {self.id}")
        if self.upload: # 3.上传
            cmd.append(f"cp -r {self.out}/{self.id}.zip {self.upload}")
            cmd.append(f"cp -r {self.p_report} {self.upload}")
        #log
        cmd.append(f"cat {self.log}/*sh.e >{self.out}/log.e.txt")
        common.cmd2shell(cmd,sh)
        return sh


    def run(self):
        """运行"""
        self.make_dir()
        sh = self.run_report()
        if not self.run_flag:
            common.prun((sh ,self.log))


#  ===== Option ======================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i','--pipeline_yaml',type=click.Path(),help="分析目录下的pipeline.yaml文件")
@click.option('--norun',is_flag=True,help="是否只生成命令，不运行程序,默认运行")
def main(pipeline_yaml,norun):
    project = Report(
            pipeline_yaml = pipeline_yaml,
            run_flag = norun
            )
    project.run()



if __name__ == "__main__":
    main()