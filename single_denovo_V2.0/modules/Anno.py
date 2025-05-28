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

class Anno(Pipe):
    def __init__(self,run_flag=False,**kwargs):
        pipeline_dic = yaml.safe_load(open(kwargs['pipeline_yaml'] ))
        self.__dict__.update(pipeline_dic)
        self.run_flag = run_flag

    def creat_env(self):
        os.makedirs(self.p_VFDB ,exist_ok=True)
        os.makedirs(self.p_drug ,exist_ok=True)
        os.makedirs(self.p_eggmapper ,exist_ok=True)
        os.makedirs(self.p_CAZy ,exist_ok=True)
        os.makedirs(self.p_swissprot ,exist_ok=True)
        os.makedirs(self.p_plasmid ,exist_ok=True)


    def run_VFDB_anno(self):
        sh  = os.path.join(self.p_sh,"5.1vfdb.sh")
        cmd=[]
        cmd.append(f"""
#VFDB
source {params.ACTIVATE} denovo
cd {self.p_VFDB}
time blastn -num_threads {int(th/2)} -evalue 1e-10 -query {self.predict_fa} -db {params.db_VFDB_nt}/VFDB -outfmt '{params.blast_opt}' -perc_identity 90 -out {self.p_VFDB}/VFDB.blast
{params.python3} {BIN}/4.anno/parse_VFDB.py -1 {self.p_VFDB}/VFDB.blast -2 {self.merge_file} --ref {params.db_VFDB}/merge_VFDB.tab

# Output fasta
awk 'NR>1' {self.p_VFDB}/Virulence.txt |cut -f1 |sort |uniq >{self.p_VFDB}/Virulence_gene.id
awk 'NR>1' {self.p_VFDB}/Virulence.txt |cut -f2 |sort |uniq >{self.p_VFDB}/Virulence_ref.id
{params.seqkit} grep -f {self.p_VFDB}/Virulence_gene.id {self.predict_fa} >{self.p_VFDB}/Virulence_gene.fa
{params.seqkit} grep -f {self.p_VFDB}/Virulence_ref.id {params.db_VFDB_nt}/VFDB_setB_nt.fas >{self.p_VFDB}/Virulence_ref.fa
""")
        common.cmd2shell(cmd,sh)
        return sh



    def run_drug_anno(self):
        sh  = os.path.join(self.p_sh,"5.2drug.sh")
        cmd =[] 
        cmd.append(f"""
# rgi
source {params.ACTIVATE} rgi
cd  {self.p_drug}
rgi main -n {int(th/3)} --input_type contig --clean -i {self.predict_fa} --output_file {self.p_drug}/drug
{params.python3} {BIN}/4.anno/parse_rgi.py -1 {self.p_drug}/drug.txt -2 {self.merge_file} --ref {params.db_CARD}/card_database_v3.2.8_all.fasta.tab

####
# resfinder
source {params.ACTIVATE} resfinder
run_resfinder.py -ifa {self.predict_fa} --db_path_res {params.db_resfinder} -s "Other"  -o {self.p_drug} -acq -l 0.8 -t 0.8
{params.python3} {BIN}/4.anno/parse_resfinder.py -i {self.p_drug}/ResFinder_results_tab.txt --db {params.db_resfinder}/phenotypes.txt
""")
        common.cmd2shell(cmd,sh)
        return sh


    def run_eggmapper_anno(self):
        sh  = os.path.join(self.p_sh,"5.3eggmapper.sh")
        cmd =[] 
        Kind = 'bact' if self.kind == 'bacteria' else 'fungi'
        cmd.append(f"""#Eggmapper
source {params.ACTIVATE} denovo
cd {self.p_eggmapper} 
time emapper.py -m diamond --cpu {int(cpu*2/3)} -d {Kind} --override -i {self.predict_faa} --output {self.id} --output_dir {self.p_eggmapper} --data_dir {params.db_eggnog}
if [ -f "{self.p_eggmapper}/{self.id}.emapper.annotations" ]; then
    {params.python3} {BIN}/4.anno/EggnogParser.py -n {self.p_eggmapper}/{self.id}.emapper.annotations -f {self.predict_faa} -o {self.p_eggmapper}
    # COG
    {params.Rscript} {BIN}/4.anno/draw_COG.R {self.p_eggmapper}/COG/all.COG.class.xls {self.p_eggmapper}/COG
    
    # GO
    python3 {BIN}/4.anno/GO_anno.py -i {self.p_eggmapper}/{self.id}.emapper.annotations -db {params.db_GO}/go-basic.obo
    {params.Rscript} {BIN}/4.anno/draw_GO.R {self.p_eggmapper}/GO/GO_anno_stats.xls {self.p_eggmapper}/GO

    # KEGG
    python3 {BIN}/4.anno/KEGG_anno.py -i {self.p_eggmapper}/{self.id}.emapper.annotations -db {params.db_KEGG}
    {params.Rscript} {BIN}/4.anno/draw_KEGG.R {self.p_eggmapper}/KEGG/KEGG_anno.stat.txt {self.p_eggmapper}/KEGG
fi
rm -rf {{emappertmp_dmdn_*,Rplots.pdf}}
""")
        common.cmd2shell(cmd,sh)
        return sh


    def run_CAZy_anno(self):
        sh  = os.path.join(self.p_sh,"5.4CAZy.sh")
        p_des = f"{params.db_CAZY}/CAZyDB.07302020.fam-activities-hmm.txt"
        cmd =[] 
        cmd.append(f"""#CAZy
source {params.ACTIVATE} denovo
cd {self.p_CAZy}
time run_dbcan {self.predict_faa} protein --db_dir {params.db_CAZY}  --out_dir {self.p_CAZy}
csvtk join -t -L -f 1 {self.p_CAZy}/hmmer.out {p_des} |csvtk round -t -f 'Coverage' -n 2 |csvtk mutate -t -f 1 -p '^(\D+)' -n 'kind' >{self.p_CAZy}/CAZY.txt
csvtk freq -t -f 'kind' {self.p_CAZy}/CAZY.txt > {self.p_CAZy}/CAZY.kind_stat.txt
{params.Rscript} {BIN}/4.anno/draw_CAZy.R {self.p_CAZy}/CAZY.kind_stat.txt {self.p_CAZy}
""")
        common.cmd2shell(cmd,sh)
        return sh


#     def run_pfam_anno(self):
#         """运行时间长，取消运行"""
#         sh  = os.path.join(self.p_sh,"4.5pfam.sh")
#         cmd =[] 
#         cmd.append(f"""#pfam
# source {params.ACTIVATE} denovo
# cd {self.p_pfam}
# time pfam_scan.pl -cpu {int(cpu*2/3)} -fasta {self.predict_faa} -dir {params.Pfam} -outfile {self.p_pfam}/pfam.out
# grep -v '^#' {self.p_pfam}/pfam.out >{self.p_pfam}/pfam.txt
# sed -i -e '1i seq_id\\talignment_start\\talignment_end\\tenvelope_start\\tenvelope_end\\thmm_acc\\thmm_name\\ttype\\thmm_start\\thmm_end\\thmm_length\\tbit_score\\tE-value\\tsignificance\\tclan' -e '2d' {self.p_pfam}/pfam.txt
# """)
#         common.cmd2shell(cmd,sh)
#         return sh

    def run_swiss_prot_anno(self):
        Kind = self.kind.capitalize()
        sh  = os.path.join(self.p_sh,"5.5swiss_prot.sh")
        cmd =[] 
        cmd.append(f"{params.python3} {BIN}/4.anno/UniProt/swissprot_db.py analysing -i {self.predict_faa} -d {params.db_swissprot} --softw_align diamond --db_select {Kind} -o {self.p_swissprot}/swissprot_result.tsv")
        common.cmd2shell(cmd,sh)
        return sh

    def run_plasmid_anno(self):
        sh  = os.path.join(self.p_sh,"5.6plasmid.sh")
        cmd =[] 
        cmd.append(f"""
source {params.ACTIVATE} denovo
time blastn -num_threads {int(th/2)} -evalue 1e-5 -max_target_seqs 5000 -query {self.predict_fa} -db {params.db_plasmid}/plsdb.fna -outfmt '{params.blast_opt}' -perc_identity 90 -out {self.p_plasmid}/plasmid.blast
python3 {BIN}/4.anno/filter_PLSDB_blast.py  -1 {params.db_plasmid}/plsdb.anno.tab -2 {self.p_plasmid}/plasmid.blast
""")
        common.cmd2shell(cmd,sh)
        return sh

    def run(self):
        self.creat_env()
        sh1 = self.run_drug_anno()
        sh2 = self.run_VFDB_anno()
        sh3 = self.run_eggmapper_anno()
        sh4 = self.run_CAZy_anno()
        sh5 = self.run_swiss_prot_anno()
        sh6 = self.run_plasmid_anno()
        sh = [sh1,sh2,sh3,sh4,sh5,sh6]
        if not self.run_flag:
            common.mul_pool(sh, self.log)



##### Option ##############################################################################
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i','--pipeline_yaml',type=click.Path(),help="分析目录下的pipeline.yaml文件")
@click.option('--norun',is_flag=True,help="是否只生成命令，不运行程序,默认运行")
def main(pipeline_yaml,norun):
    project = Anno(
            pipeline_yaml = pipeline_yaml,
            run_flag = norun
            )
    project.run()



if __name__ == "__main__":
    main()