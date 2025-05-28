#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/14 9:05
# @Last Modified by:   Ming
# @Last Modified time: 2022/7/14 15:22
import logging
import multiprocessing
import os
from pathlib import Path
from subprocess import run

import yaml

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
# Some Setting
PATH = os.path.dirname(os.path.abspath(__file__))
BINPATH = os.path.realpath(os.path.join(PATH, "../bin"))
CONFIGPATH = os.path.join(PATH, "../config")
f_software = os.path.join(CONFIGPATH, "software.yml")
f_env = os.path.join(CONFIGPATH, "env.yml")
f_database = os.path.join(CONFIGPATH, "database.yml")


class PTyping(object):
    """
    分型流程

    TODO: 逻辑整理
    """

    def __init__(self, fasta, database, out):
        """
        Init the object

        :param fasta: The input fasta genome file
        :param database: The database to use
        :param out: The output dir
        """
        self.env = yaml.safe_load(open(f_env, 'r').read())
        self.software = yaml.safe_load(open(f_software, 'r').read())

        self.fasta = os.path.abspath(fasta)
        self.db_name = database
        self.database = yaml.safe_load(open(f_database, 'r').read())["typing"]
        self.out = os.path.abspath(out)
        os.makedirs(self.out, exist_ok=True)
        self.script = os.path.join(self.out, "typing.sh")
        self.log = os.path.join(self.out, "typing.log")

        # 机器信息
        self.cpu = multiprocessing.cpu_count()

    def run(self):
        """

        """
        sub_run = {"沙门氏菌属": self._salmonella(),
                   "大肠埃希菌": self._ecoli(),
                   "人类免疫缺陷病毒": self._hiv(),
                   "炭疽杆菌属": self._bacillus(),
                   "腮腺炎病毒": self._genetree("Mumps orthorubulavirus", f"{self.database}/mor/mor.aln"),
                   "风疹病毒": self._genetree("Rubella virus", f"{self.database}/rvi/rvi.aln"),
                   "麻疹病毒": self._genetree("Measles morbillivirus", f"{self.database}/mmo/mmo.aln"),
                   "西尼罗病毒": self._genetree("West Nile virus", f"{self.database}/wnv/wnv.aln"),
                   "蜱传脑炎病毒": self._genetree("Tick-borne encephalitis virus", f"{self.database}/tev/tev.aln"),
                   "札如病毒": self._genetree("Sapporo virus", f"{self.database}/svi/svi.aln"),
                   "人偏肺病毒": self._genetree("hMPV", f"{self.database}/hmpv/hmpv.aln"),
                   "裂谷热病毒": self._genetree("Rift Valley fever virus", f"{self.database}/rvfv/rvfv.aln"),
                   "肺炎衣原体": self._genetree("Chlamydia pneumoniae", f"{self.database}/cpn/cpn.aln"),
                   "土拉弗朗西斯菌": self._ani("Francisella tularensis", f"{self.database}/ftu"),
                   "诺如病毒": self._noro(),
                   "登革热病毒": self._genetree_ql("Dengue virus", f"{self.database}/denv/denv.gb"),
                   "汉坦病毒属": self._genetree_ql("orthohantavirus", f"{self.database}/hanta/hanta.gb"),
                   "肠道病毒通用型": self._genetree_ql("enterovirus", f"{self.database}/entv/entv.gb"),
                   "人腺病毒": self._genetree_ql("adenovirus", f"{self.database}/adenovirus/adenovirus.gb"),
                   "呼吸道合胞病毒": self._genetree_ql("Respiratory Syncytial Virus", f"{self.database}/rsv/rsv.gb"),
                   "轮状病毒": self._genetree_ql("rotavirus", f"{self.database}/rotav/rotav.gb"),
                   "乙型脑炎病毒": self._genetree_ql("Japanese encephalitis virus", f"{self.database}/jev/jev.gb"),
                   "埃博拉病毒": self._genetree_ql("Ebolavirus", f"{self.database}/ebola/ebola.gb"),
                   "狂犬病病毒": self._genetree_ql("Lyssavirus rabies", f"{self.database}/rabv/rabv.gb"),
                   "寨卡病毒": self._genetree_ql("zika virus", f"{self.database}/zika/zika.gb"),
                   "人副流感病毒": self._genetree_ql("Orthorubulavirus", f"{self.database}/hpiv/hpiv.gb"),
                   "乙肝病毒": self._genetree_ql("Hepatitis B Virus", f"{self.database}/hbv/hbv.gb"),
                   "丙肝病毒": self._genetree_ql("Hepatitis C Virus", f"{self.database}/hcv/hcv.gb"),
                   "人疱疹病毒": self._genetree_ql("Herpesviridae", f"{self.database}/hhv/hhv.gb"),
                   "人博卡病毒": self._genetree_ql("Human Bocavirus", f"{self.database}/hbov/hbov.gb"),
                   "人星状病毒": self._genetree_ql("Mamastrovirus", f"{self.database}/hastv/hastv.gb"),
                   "布尼亚病毒属": self._genetree_ql("orthobunyavirus", f"{self.database}/bunya/bunya.gb"),
                   "人鼻病毒": self._genetree_ql("human rhinovirus", f"{self.database}/hrv/hrv.gb"),
                   "基孔肯雅病毒": self._genetree_ql("Chikungunya virus", f"{self.database}/chikv/chikv.gb"),
                   "尼帕病毒": self._genetree_ql("Nipah henipavirus", f"{self.database}/niv/niv.gb"),
                   "新疆出血热病毒": self._genetree_ql("Crimean-Congo hemorrhagic fever virus",
                                                       f"{self.database}/ccfhv/ccfhv.gb"),
                   "拉沙热病毒": self._genetree_ql("lassa virus", f"{self.database}/lassa/lassa.gb"),
                   "亨德拉病毒": self._genetree_ql("Hendra virus", f"{self.database}/hev/hev.gb"),
                   "猴痘病毒": self._mpxv(),
                   "黄热病毒": self._genetree_ql("Yellow Fever virus", f"{self.database}/yfv/yfv.gb"),
                   "甲肝病毒": self._genetree_ql("Hepatitis A Virus", f"{self.database}/hav/hav.gb"),
                   "甲型流感病毒": self._flu(),
                   "乙型流感病毒": self._fluB(),
                   "丙型流感病毒": self._flu(),
                   "唐菖蒲伯克霍尔德氏菌": self._ani("Burkholderia gladioli", f"{self.database}/bgl"),
                   "脑膜炎奈瑟菌": self._neisseria(),
                   "霍乱弧菌": self._vibrio_cholerae(),
                   "志贺氏菌": self._shigella(),
                   "单核细胞增生李斯特菌": self._listeria("listeria"),
                   "流感嗜血杆菌": self._haemophilus_influenzae(),
                   "肺炎克雷伯杆菌": self._kpn(),
                   "副溶血弧菌": self._vibrio_parahaemolyticus()}

        self.command = sub_run[self.db_name]
        with open(self.script, 'w') as OUT:
            print(self.command, file=OUT)
        run(f"{self.software['shell']} {self.script} 2>&1 > {self.log}",
            shell=True)

    def to_upload(self, d_upload):
        """
        link the result file to dir upload
        """
        cmd = f"cp -rf {self.out}/result.tsv {d_upload}/typing.tsv"
        run(cmd, shell=True)

    @property
    def result(self):
        """
        获取血清型分型结果
        """
        self.f_result = Path(f"{self.out}/result.tsv").absolute()
        if self.f_result.exists():
            sub_run = {"沙门氏菌属": self._salmonella_result(),
                       "大肠埃希菌": self._ecoli_result(),
                       "人类免疫缺陷病毒": self._hiv_result(),
                       "炭疽杆菌属": self._bacillus_result(),
                       "腮腺炎病毒": self._genetree_result(),
                       "风疹病毒": self._genetree_result(),
                       "麻疹病毒": self._genetree_result(),
                       "西尼罗病毒": self._genetree_result(),
                       "蜱传脑炎病毒": self._genetree_result(),
                       "札如病毒": self._genetree_result(),
                       "人偏肺病毒": self._genetree_result(),
                       "裂谷热病毒": self._genetree_result(),
                       "肺炎衣原体": self._genetree_result(),
                       "土拉弗朗西斯菌": self._ani_result(),
                       "诺如病毒": self._noro_result(),
                       "登革热病毒": self._genetree_ql_result(),
                       "汉坦病毒属": self._genetree_ql_result(),
                       "肠道病毒通用型": self._genetree_ql_result(),
                       "人腺病毒": self._genetree_ql_result(),
                       "呼吸道合胞病毒": self._genetree_ql_result(),
                       "轮状病毒": self._genetree_ql_result(),
                       "乙型脑炎病毒": self._genetree_ql_result(),
                       "埃博拉病毒": self._genetree_ql_result(),
                       "狂犬病病毒": self._genetree_ql_result(),
                       "寨卡病毒": self._genetree_ql_result(),
                       "人副流感病毒": self._genetree_ql_result(),
                       "乙肝病毒": self._genetree_ql_result(),
                       "丙肝病毒": self._genetree_ql_result(),
                       "人疱疹病毒": self._genetree_ql_result(),
                       "人博卡病毒": self._genetree_ql_result(),
                       "人星状病毒": self._genetree_ql_result(),
                       "布尼亚病毒属": self._genetree_ql_result(),
                       "人鼻病毒": self._genetree_ql_result(),
                       "基孔肯雅病毒": self._genetree_ql_result(),
                       "尼帕病毒": self._genetree_ql_result(),
                       "新疆出血热病毒": self._genetree_ql_result(),
                       "拉沙热病毒": self._genetree_ql_result(),
                       "亨德拉病毒": self._genetree_ql_result(),
                       "猴痘病毒": self._mpxv_result(),
                       "黄热病毒": self._genetree_ql_result(),
                       "甲肝病毒": self._genetree_ql_result(),
                       "甲型流感病毒": self._flu_result(),
                       "乙型流感病毒": self._fluB_result(),
                       "丙型流感病毒": self._flu_result(),
                       "唐菖蒲伯克霍尔德氏菌": self._ani_result(),
                       "脑膜炎奈瑟菌": self._neisseria_result(),
                       "霍乱弧菌": self._vibrio_cholerae_result(),
                       "志贺氏菌": self._shigella_result(),
                       "单核细胞增生李斯特菌": self._listeria_result(),
                       "流感嗜血杆菌": self._haemophilus_influenzae_result(),
                       "肺炎克雷伯杆菌": self._kpn_result(),
                       "副溶血弧菌": self._vibrio_parahaemolyticus_result()}
            return sub_run[self.db_name]
        else:
            return None

    def _flu(self):
        """
        流感病毒
        """
        LABEL = self.software['label']
        cmd = f"""# 流感
set -e
cd {self.out}
bash {LABEL} {self.fasta} result irma-FLU-v2
sed -e '1s,VIRUS STRAIN,序列名称,g' -e '1s,CLADE,分型预测,g' result_final.txt > result.tsv
"""
        return cmd

    def _flu_result(self):
        """"""
        h = "H?"
        n = "N?"
        with open(self.f_result, 'r') as IN:
            for line in IN:
                arr = line.strip().split("\t")
                if arr[1].startswith("A_HA"):
                    h = arr[1].strip().split('_')[-1]
                if arr[1].startswith("A_NA"):
                    n = arr[1].strip().split('_')[-1]
        return f"{h}{n}"

    def _fluB(self):
        """
        乙型流感病毒
        """
        LABEL = self.software['label']
        cmd = f"""# 乙型流感病毒
set -e
cd {self.out}
bash {LABEL} {self.fasta} all irma-FLU-v2
for i in B_HAv2019 B_MPv2016 B_NAv2016 B_PAv2016 B_NSv2016 B_PAv2016 B_PB1v2016 B_PB2v2016
do
  bash {LABEL} {self.fasta} $i $i
done
python {BINPATH}/merge_fluB.py --flu all_final.txt --segments B_HAv2019_final.txt,B_MPv2016_final.txt,B_NAv2016_final.txt,B_PAv2016_final.txt,B_NSv2016_final.txt,B_PAv2016_final.txt,B_PB1v2016_final.txt,B_PB2v2016_final.txt --out result.tsv
"""
        return cmd

    def _fluB_result(self):
        """

        """
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split("\t")
                return arr[-1]

    def _salmonella(self):
        """
        沙门氏菌
        """
        cmd = f"""# Salmonella
set -e
source {self.software["activate"]} {self.env["typetrace"]}
sistr -f tab -o {self.out}/tmp {self.fasta}
awk -F'\t' -vOFS="\t" 'BEGIN{{print "亚种", "H1", "H2", "O_antigen", "血清群", "血清型"}}NR>1{{print $6,$9,$10,$11,$12,$13}}' {self.out}/tmp.tab > {self.out}/result.tsv
"""
        return cmd

    def _salmonella_result(self):
        """

        """
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split("\t")
                return arr[-1]

    def _ecoli(self):
        """
        大肠埃希菌

        TODO: 重写sero_ecoli.py
        """
        path_vfdb = yaml.safe_load(open(f_database, 'r').read())["vfdb"]
        res = f"""# Escherichia coli
set -e
source {self.software["activate"]} {self.env["typetrace"]}
python {BINPATH}/sero_ecoli.py --extented_output -i {self.fasta} -p {self.database}/ecoli -o {self.out}
blastn -db {path_vfdb}/VFDB_setB_nt/VFDB -query {self.fasta} -num_threads 10 -outfmt '6 qseqid stitle pident' > {self.out}/vfdb_blastn_res.tab
{self.software['rscript']} --vanilla {BINPATH}/find_subtype_marker_gene_ecoli_pathotype.R {self.out}/vfdb_blastn_res.tab 95 {self.out}/results_tab.tsv
ln -sf {self.out}/results_tab.tsv {self.out}/result.tsv 
"""
        return res

    def _ecoli_result(self):
        """

        """
        h = "H?"
        o = "O?"
        set_o = set()
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split("\t")
                if line.startswith("H_type"):
                    h = arr[2]
                if line.startswith("O_type"):
                    set_o.add(arr[2])
                if len(set_o) > 0:
                    o = '|'.join(list(set_o))

        return f"{o}:{h}"

    def _hiv(self):
        """
        人类免疫缺陷病毒
        """
        res = f"""# Human immunodeficiency virus
set -e
{self.software['rscript']} {BINPATH}/genome2segment.R {self.fasta} {self.out}/segment.fasta 300 100
{self.software['blastn']} -query {self.out}/segment.fasta -db {self.database}/hiv/hiv -out {self.out}/segment.blast.res -outfmt 6
if [[ $(cat {self.out}/segment.blast.res) == '' ]];
then
    echo -e "分型\t亚型\nNA\tNA" > {self.out}/result.tsv
else
    awk -F '\\t' '$3>30' {self.out}/segment.blast.res|awk -F'_' -vOFS="\\t" '{{print ($5-1)/100+1,$0}}'| sort -u -k1n -k13nr > {self.out}/sorted_segment.tsv
    awk '{{if ($2>a[$1]){{a[$1]=$2;print}}}}' {self.out}/sorted_segment.tsv | sort -n > {self.out}/max_score_segment.tsv
    max=$(less -S {self.out}/sorted_segment.tsv | sort -k3 | awk -F '\\t' '{{a[$3] += $11}}END{{for (i in a)print i, a[i]}}' | sort -nr -k2 | head -1 | cut -d ' ' -f1)
    echo $max | awk -F'.' '{{print $1}}' | xargs -I {{}} grep {{}} {self.database}/hiv/metadata.tsv | cut -f 2- | sed 1i"分型\t亚型" > {self.out}/result.tsv
fi  
"""
        return res

    def _hiv_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split("\t")
                return arr[-1]

    def _bacillus(self):
        """
        炭疽杆菌
        """
        res = f"""# Bacillus anthraci
set -e
source {self.software["activate"]} {self.env["typetrace"]}
mkdir -p {self.out}/before_QC
{self.software['wgsim']} -1 200 -2 200 -r 0 -R 0 -X 0 -e 0 {self.fasta} {self.out}/before_QC/query_1.fastq {self.out}/before_QC/query_2.fastq
mkdir -p {self.out}/after_QC
fastp -i {self.out}/before_QC/query_1.fastq -I {self.out}/before_QC/query_2.fastq -o {self.out}/after_QC/query_1.fastq -O {self.out}/after_QC/query_2.fastq -j {self.out}/after_QC/query.json -h {self.out}/after_QC/query.html
mkdir -p {self.out}/bwa
{self.software['bwa']} mem -t 4 {self.database}/bacillus/bacillus {self.out}/after_QC/query_1.fastq {self.out}/after_QC/query_2.fastq > {self.out}/bwa/query.sam
{self.software['samtools']} sort -o {self.out}/bwa/query.bam {self.out}/bwa/query.sam
{self.software['samtools']} index {self.out}/bwa/query.bam
gatk MarkDuplicates --QUIET true --VERBOSITY ERROR -I {self.out}/bwa/query.bam -O {self.out}/bwa/query.rmdup.bam -M {self.out}/bwa/query_dup_metrics.txt
gatk AddOrReplaceReadGroups --QUIET true --VERBOSITY ERROR -I {self.out}/bwa/query.rmdup.bam -O {self.out}/bwa/query.addgp.bam -SM query -LB query -PL ILLUMINA -PU unit1
{self.software['samtools']} index {self.out}/bwa/query.addgp.bam
mkdir -p {self.out}/variant
gatk HaplotypeCaller --QUIET true --verbosity ERROR -R {self.database}/bacillus/bacillus.fasta -I {self.out}/bwa/query.addgp.bam -O {self.out}/variant/query.g.vcf.gz -ERC GVCF -ploidy 1
mkdir -p {self.out}/combined_variant
gatk CombineGVCFs --QUIET true --verbosity ERROR -R {self.database}/bacillus/bacillus.fasta -V {self.out}/variant/query.g.vcf.gz -V {self.database}/bacillus/combine.g.vcf.gz -O {self.out}/combined_variant/query.g.vcf.gz
gatk GenotypeGVCFs --QUIET true --verbosity ERROR -R {self.database}/bacillus/bacillus.fasta -V {self.out}/combined_variant/query.g.vcf.gz -O {self.out}/variant/query.vcf.gz
vcftools --gzvcf {self.out}/variant/query.vcf.gz --remove-indels --recode --recode-INFO-all --out {self.out}/variant/query.SNP
mkdir -p {self.out}/phylip
python {BINPATH}/vcf2phylip.py -i {self.out}/variant/query.SNP.recode.vcf --output-prefix {self.out}/phylip/query -m 13
mkdir -p {self.out}/tree
{self.software['fasttree']} -nt {self.out}/phylip/query.min13.phy > {self.out}/tree/query.tree
mkdir -p {self.out}/final_results
{self.software['rscript']} --vanilla {BINPATH}/plottree_subtyping.R {self.out}/tree/query.tree {self.out} query
"""
        return res

    def _bacillus_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split("\t")
                return arr[-1]

    def _genetree(self, name, ref):
        """
        叶辉这边利用基因建树的通用流程

        TODO: 速度太慢，需要寻找新的方法
        """
        d_database = os.path.dirname(ref)
        cmd = f"""# GeneTree for {name}
set -e 
source {self.software["activate"]} {self.env["typetrace"]}
{self.software['python']} {BINPATH}/seq_merge.py --fasta {self.fasta} --out {self.out}/query.fasta
# awk '{{if (/^>/) print $1;else print $0}}' \"{self.fasta}\" > {self.out}/query.fasta
mafft --quiet --add {self.out}/query.fasta --reorder {ref} > {self.out}/merge.aln
bmge -i {self.out}/merge.aln -t DNA -m DNAPAM100:3 -b 1 -o {self.out}/merge.mafa
iqtree --quiet -keep-ident -s {self.out}/merge.mafa --redo -m MFP -bb 1000 -nt {int(self.cpu / 2)} --prefix {self.out}/tree
id=`grep ">" {self.out}/query.fasta | awk -F" " '{{print $1}}'| sed s'/\./_/'|sed s'/>//g'`
perl {BINPATH}/genetree.pl {self.out}/tree.mldist $id {d_database}/metadata.tsv {self.out}
"""
        return cmd

    def _genetree_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _genetree_ql(self, name, ref):
        """
        千里这边利用基因建树的通用流程
        """
        cmd = f"""# GeneTree for {name}
set -e 
source {self.software["activate"]} {self.env["typetrace"]}
cd {self.out}
python {BINPATH}/subtyping_virus.py -p subtype -i {self.fasta} -p subtype -db {ref} -c Serotype -b 100 -t 0.1
"""
        return cmd

    def _genetree_ql_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _ani(self, name, ref):
        """
        ANI 分型方法
        """
        d_database = os.path.dirname(ref)
        cmd = f"""# ANI for {name}
set -e 
source {self.software["activate"]} {self.env["typetrace"]}
cp -rf {ref}/*.fna {self.out}
cp {self.fasta} {self.out}/query.fna
average_nucleotide_identity.py -i {self.out} -o {self.out}/ANI -m ANIm -f
{self.software['rscript']} --vanilla {BINPATH}/ani_subtyping.R {self.out}/ANI/ANIm_percentage_identity.tab query {self.out}/result.tsv 
"""
        return cmd

    def _ani_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _neisseria(self):
        """
        脑膜炎奈瑟菌
        """
        cmd = f"""# Neisseria meningitidis
set -e
mkdir -p {self.out}/genome
cp {self.fasta} {self.out}/genome/
python {BINPATH}/characterize_neisseria_capsule.py -d {self.out}/genome/ --database {self.database}/neisseria -t {self.cpu // 2} -o {self.out}
cut -f 1,2 {self.out}/serogroup/serogroup_predictions.tab > {self.out}/result.tsv
"""
        return cmd

    def _neisseria_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _noro(self):
        """
        诺如病毒
        """
        cmd = f"""# Norovirus
set -e
source {self.software["activate"]} {self.env["typetrace"]}
cd {self.out}
python {BINPATH}/noro.py -p subtype -db {self.database}/noro/noro.gb  -b 100 -t 0.1 -m {self.database}/noro/metadata.tsv -a -g ""RdRp#dependent#polymerase"" ""VP1#Capsid""  -c P_Type -i {self.fasta}
"""
        return cmd

    def _noro_result(self):
        result = set()
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split("\t")
                if line.startswith("RdRp"):
                    result.add(arr[1].strip().split('.')[0])
                if line.startswith("VP1"):
                    result.add(arr[1].strip().split('.')[0])
        return '|'.join(list(result))

    def _vibrio_cholerae(self):
        """
        霍乱弧菌
        """
        cmd = f"""# vibrio cholerae
set -e
mkdir -p {self.out}/simulated_reads
{self.software['wgsim']} -1 200 -2 200 -r 0 -R 0 -X 0 -e 0 {self.fasta} {self.out}/simulated_reads/1.fastq {self.out}/simulated_reads/2.fastq
mkdir -p {self.out}/alignment_res
{self.software['bowtie2']} -x {self.database}/Vcholerae/wbeO1 -1 {self.out}/simulated_reads/1.fastq -2 {self.out}/simulated_reads/2.fastq > {self.out}/alignment_res/wbeO1.sam
{self.software['bowtie2']} -x {self.database}/Vcholerae/wbfO139 -1 {self.out}/simulated_reads/1.fastq -2 {self.out}/simulated_reads/2.fastq > {self.out}/alignment_res/wbeO139.sam
{self.software['bowtie2']} -x {self.database}/Vcholerae/ctx -1 {self.out}/simulated_reads/1.fastq -2 {self.out}/simulated_reads/2.fastq > {self.out}/alignment_res/ctx.sam
{self.software['samtools']} sort -o {self.out}/alignment_res/wbeO1.bam {self.out}/alignment_res/wbeO1.sam
{self.software['samtools']} sort -o {self.out}/alignment_res/wbeO139.bam {self.out}/alignment_res/wbeO139.sam
{self.software['samtools']} sort -o {self.out}/alignment_res/ctx.bam {self.out}/alignment_res/ctx.sam
{self.software['samtools']} index {self.out}/alignment_res/wbeO1.bam
{self.software['samtools']} index {self.out}/alignment_res/wbeO139.bam
{self.software['samtools']} index {self.out}/alignment_res/ctx.bam
{self.software['samtools']} coverage {self.out}/alignment_res/wbeO1.bam > {self.out}/alignment_res/wbeO1.coverage
{self.software['samtools']} coverage {self.out}/alignment_res/wbeO139.bam > {self.out}/alignment_res/wbeO139.coverage
{self.software['samtools']} coverage {self.out}/alignment_res/ctx.bam > {self.out}/alignment_res/ctx.coverage
{self.software['rscript']} --vanilla {BINPATH}/Vibrio_cholerae_serogroup_find_coverage.R {self.fasta} {self.out}
{self.software['rscript']} --vanilla {BINPATH}/Vibrio_cholerae_ctx_find_coverage.R {self.fasta} {self.out}
"""
        return cmd

    def _vibrio_cholerae_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _shigella(self):
        """
        志贺氏菌
        """
        cmd = f"""# Shigella
set -e
source {self.software["activate"]} {self.env["shigatpyer"]}
mkdir -p {self.out}/simulated_reads
{self.software['wgsim']} -1 200 -2 200 -r 0 -R 0 -X 0 -e 0 {self.fasta} {self.out}/simulated_reads/1.fastq {self.out}/simulated_reads/2.fastq
cd {self.out}
shigatyper --R1 {self.out}/simulated_reads/1.fastq --R2 {self.out}/simulated_reads/2.fastq -n shiga
ln -sf {self.out}/shiga.tsv {self.out}/result.tsv
"""
        return cmd

    def _shigella_result(self):
        """"""
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _listeria(self, typing_scheme):
        """
        单核细胞增生李斯特菌
        """
        blastdb = os.path.join(self.database, "listeria/blast/mlst.fa")
        datadir = os.path.join(self.database, "listeria/pubmlst")
        cmd = f"""# Listeria monocytogenes 
set -e
source {self.software["activate"]} {self.env["typetrace"]}
mlst --nopath --legacy --quiet --scheme {typing_scheme} -blastdb {blastdb} -datadir {datadir} \"{self.fasta}\" > {self.out}/mlst_listeria_serotype.tmp
unknown_types=`less {self.out}/mlst_listeria_serotype.tmp | sed 1d | cut -d $'\\t' -f4`
if [[ $unknown_types =~ '-' ]]
then
    sero_id_cols='Sero_id'
    sero_type_cols='Sero_type'
    for i in `echo unknown_types`
    do
        result=`less {datadir}/{typing_scheme}/{typing_scheme}.txt`
        gene_types=`less {self.out}/mlst_listeria_serotype.tmp | cut -d $'\\t' -f5,6,7,8,9|sed 1d`
        ncol=1
        for j in `echo $gene_types`
        do
            ncol=`expr $ncol + 1`
            [ $j == '-' ] && continue
            result=`echo -e \"$result\" | awk -F '\\t' -v x=$j '$'$ncol'==x'`
        done
        n=`echo -e \"$result\" | wc -l`
        if [ $n -eq 1 ]
        then
            sero_id=`echo $result|cut -d ' ' -f1` && sero_type=`echo $result|cut -d ' ' -f7`;else sero_id='-' && sero_type='-'
        fi
        sero_id_cols=$sero_id_cols'\\n'$sero_id
        sero_type_cols=$sero_type_cols'\\n'$sero_type
    done
    paste -d $'\\t' <(cut -d $'\\t' -f1,2,3 {self.out}/mlst_listeria_serotype.tmp) <(echo -e \"$sero_id_cols\") <(cut -d $'\\t' -f5,6,7,8,9 {self.out}/mlst_listeria_serotype.tmp) <(echo -e \"$sero_type_cols\") > {self.out}/result.tmp
fi
awk -F'\\t' -vOFS="\\t" 'BEGIN{{print "序列名称", "血清型"}}NR==2{{print "Query", $10}}' {self.out}/result.tmp > {self.out}/result.tsv
"""
        return cmd

    def _listeria_result(self):
        """

        """
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _haemophilus_influenzae(self):
        """
        流感嗜血杆菌
        """
        cmd = f"""# Haemophilus influenzae
set -e
source {self.software["activate"]} {self.env["hicap"]}
ln -sf \'{self.fasta}\' {self.out}/query.fasta
hicap -d {self.database}/hin -q {self.out}/query.fasta -o {self.out}
ln -sf {self.out}/query.tsv {self.out}/result.tsv
"""
        return cmd

    def _haemophilus_influenzae_result(self):
        """

        """
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _kpn(self):
        """
        肺炎克雷伯

        TODO: 添加软件输出结果的解释
        """
        cmd = f"""# Klebsiella pneumoniae
set -e
source {self.software["activate"]} {self.env["kleborate"]}
ln -sf \'{self.fasta}\' {self.out}/query.fasta
kleborate --resistance -o {self.out}/result.tsv -a {self.out}/query.fasta
"""
        return cmd

    def _kpn_result(self):
        """"""
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _mpxv(self):
        """
        猴痘流程
        """
        cmd = f"""# Monkeypox
set -e
source {self.software["activate"]} {self.env["nextclade"]}
ln -sf \'{self.fasta}\' {self.out}/query.fasta
nextclade run -D {self.database}/mpxv -t {self.out}/nextclade.tsv {self.out}/query.fasta
cut -f 1,2,4 {self.out}/nextclade.tsv > {self.out}/result.tsv
"""
        return cmd

    def _mpxv_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]

    def _vibrio_parahaemolyticus(self):
        """
        副溶血弧菌
        """
        cmd = f"""# Vibrio parahaemolyticus
set -e
cd {self.out}
ln -sf \'{self.fasta}\' {self.out}/query.fasta
{self.software['python']} {BINPATH}/Vibrio_parahaemolyticus_sero.py \
-f {self.out}/query.fasta \
-d {self.database}/Vibrio_parahaemolyticus/reference \
-l {self.database}/Vibrio_parahaemolyticus/length.tsv \
-o {self.out} \
--blastn {self.software['blastn']}
"""
        return cmd

    def _vibrio_parahaemolyticus_result(self):
        with open(self.f_result, 'r') as IN:
            next(IN)
            for line in IN:
                arr = line.strip().split('\t')
                return arr[1]
