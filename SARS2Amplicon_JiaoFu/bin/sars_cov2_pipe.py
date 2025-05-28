#!/usr/bin/env python

import os
import sys
import re
import argparse
import logging
from subprocess import run
import yaml
from pathlib import Path
import pdb
# from ..lib.Modules import Alignment
sys.path.append(f"{Path(__file__).resolve().parents[1]}")
from lib.Modules import Alignment


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)


def get_args():
    parser = argparse.ArgumentParser()
    parser.description="二代测序新冠扩增子分析流程, 支持单/双端FASTQ, 可选上传引物BED去引物"
    parser.add_argument("-i", "--fastq1", required=True, help="单端测序输入的FASTQ, 或双端测序的read1.")
    parser.add_argument("-I", "--fastq2", help="双端测序需要填, 双端测序的read.")
    parser.add_argument("-o", "--outdir", default="./analysis", help="结果输出的文件夹.")
    parser.add_argument("-p", "--prefix", default="sample", help="样本名称.")
    parser.add_argument("-t", "--bed", help="扩增子的靶向区域, BED文件包含参考名和位置信息.")
    parser.add_argument("--global_sars2", dest="gsars2", action="store_true", help="是否使用新冠全球库进行分析, 默认不运行.")
    parser.add_argument('--trim_software', default='ivar', choices=["ivar","fgbio","ktrim"],
                        help='去引物程序选择,与-t/--bed参数一起使用. ["ivar","fgbio","ktrim"]')
    parser.add_argument("-d", "--dryrun", action="store_true", help="是否运行程序, 如选择生成shell脚本不运行.")
    return parser.parse_args()

class Amplicon():
    """
    单样本流程
        args :: 命令行参数 Namespace
    """
    def __init__(self, args):
        self.dir_bin = os.path.dirname(__file__)
        self.dir_conf = Path(__file__).parents[1].joinpath("conf").resolve()
        self.dir_etc = Path(__file__).parents[1].joinpath("etc").resolve()
        self.assign_base_info(args)
        self.assign_software()
        self.assign_parameters()
        self.assign_databases()
        self.cmds = list()

    def assign_base_info(self, args):
        logging.info("读取输入信息")
        self.fastq1 = args.fastq1
        self.outdir = args.outdir
        self.prefix = args.prefix
        self.fastq2 = args.fastq2
        self.bed    = args.bed
        self.gsars2 = args.gsars2 # 布尔值,新冠全球库做不做
        self.trim_software = args.trim_software # 去引物软件放在参数里
        self.dryrun = args.dryrun

    def assign_software(self):
        logging.info("读取软件配置文件")
        with open(os.path.join(self.dir_conf, "software.yaml"), "rt") as fh:
            self.dict_soft = yaml.safe_load(fh)
        self.python = self.dict_soft["python"]
        self.fastqc = self.dict_soft["fastqc"]
        self.fastp = self.dict_soft["fastp"]
        self.bwa = self.dict_soft["bwa"]
        self.samtools = self.dict_soft["samtools"]
        self.Rscript = self.dict_soft["Rscript"]
        self.bedtools = self.dict_soft["bedtools"]
        self.ivar = self.dict_soft["ivar"]
        self.freebayes = self.dict_soft["freebayes"]
        self.bgzip = self.dict_soft["bgzip"]
        self.bcftools = self.dict_soft["bcftools"]
        self.snpEff = self.dict_soft["snpEff"]
        self.vcftools = self.dict_soft["vcftools"]
        self.mafft = self.dict_soft["mafft"]
        self.fasttree = self.dict_soft["fasttree"]
        self.magick = self.dict_soft["magick"]
        self.gatk = self.dict_soft["gatk"]
        self.fgbio = self.dict_soft["fgbio"]
        self.ktrim = self.dict_soft["ktrim"]

    def assign_parameters(self):
        logging.info("读取软件参数配置文件")
        with open(os.path.join(self.dir_conf, "parameters.yaml"), "rt") as fh:
            self.dict_params = yaml.safe_load(fh)
        self.thread = self.dict_params["thread"]
        self.parallel_number = self.dict_params["parallel_number"]
        self.min_cov = self.dict_params["min_cov"]
        self.min_alt_count = self.dict_params["min_alt_count"]
        self.min_alt_frac = self.dict_params["min_alt_frac"]
        self.options_fastp = self.dict_params["options_fastp"]
        # self.trim_software = self.dict_params["trim_software"]

    def assign_databases(self):
        logging.info("读取软件参数配置文件")
        with open(os.path.join(self.dir_conf, "databases.yaml"), "rt") as fh:
            self.dict_db = yaml.safe_load(fh)
        self.bwaidx = self.dict_db["bwadb_prefix"]

    def qc(self):
        """
        1. 原始数据fastqc
        2. fastp过滤
        3. fastp结果统计表格
        4. 过滤后数据fastqc [221028]
        """
        logging.info("原始数据QC")
        if self.fastq2:
            fastq2 = self.fastq2
            fastp2I = f"-I {self.fastq2}"
            fastp2O = f"-O {self.dict_dirs['qc']}/{self.prefix}.clean.2.fq"
            parse_fastp_mode = "PE"
        else:
            fastq2,fastp2I,fastp2O = "", "", ""
            parse_fastp_mode = "SE"
        self.cmds.append(f"""
# 原始数据QC
mkdir {self.dict_dirs['qc']}/before {self.dict_dirs['qc']}/after
{self.fastqc} -t {self.thread} --extract -o {self.dict_dirs['qc']}/before {self.fastq1} {fastq2} &
{self.fastp} --thread {self.thread} {self.options_fastp} \
    --adapter_fasta {self.dir_etc}/adapter.fa \
    -j {self.dict_dirs['qc']}/{self.prefix}.json -h {self.dict_dirs['qc']}/{self.prefix}.html \
    -i {self.fastq1} {fastp2I} \
    -o {self.dict_dirs['qc']}/{self.prefix}.clean.1.fq {fastp2O}
{self.python} {self.dir_bin}/fastp_json2table.py \
    {self.dict_dirs['qc']}/{self.prefix}.json \
    {self.dict_dirs['qc']}/{self.prefix}.basic.stat.txt
# {self.python} {self.dir_bin}/parse_fastp.json.py \
#     -m {parse_fastp_mode} -id {self.prefix} \
#     -i {self.dict_dirs['qc']}/{self.prefix}.json
{self.fastqc} -t {self.thread} --extract \
    {self.dict_dirs['qc']}/{self.prefix}.clean.1.fq {fastp2O.replace('-O ','')} \
    -o {self.dict_dirs['qc']}/after 
        """)

    def alignment(self):
        logging.info("比对")
        _cml = Alignment(self.__dict__).pipe()
        self.cmds.append(_cml)
        # 深度统计 
        self.cmds.append(f"""
#~ 深度统计
{self.samtools} depth -a {self.dict_dirs['map']}/{self.prefix}.bam > {self.dict_dirs['map']}/{self.prefix}.bam.depth
{self.samtools} stat {self.dict_dirs['map']}/{self.prefix}.bam \
    | grep ^SN | cut -f 2- \
    > {self.dict_dirs['map']}/{self.prefix}.bam.summary_numbers
{self.Rscript} {self.dir_bin}/genome_coverage.R \
    {self.dict_dirs['map']}/{self.prefix}.bam.depth \
    {self.dict_dirs['map']}/{self.prefix}.bam.summary_numbers \
    {self.dict_dirs['map']}/{self.prefix} &
#~ 低深度区域
{self.bedtools} genomecov -bga -ibam {self.dict_dirs['map']}/{self.prefix}.bam \
    | awk '$4<{self.min_cov}' | {self.bedtools} merge -i - \
    | awk '$3-$2>20' 1> {self.dict_dirs['map']}/{self.prefix}.lowcovmask.bed
        """)

    def freebayes_variants(self):
        self.cmds.append(f"""
{self.freebayes} -F {self.min_alt_frac} -C {self.min_alt_count} \
    -f {self.dict_dirs['consensus']}/{self.prefix}.masked.fa \
    {self.dict_dirs['map']}/{self.prefix}.bam \
    > {self.dict_dirs['variant']}/{self.prefix}.raw.vcf
        """)

    def gatk_variant_pipe(self):
        self.cmds.append(f"{self.gatk} CreateSequenceDictionary "
                        f"-R {self.dict_dirs['consensus']}/{self.prefix}.masked.fa "
                        f"-O {self.dict_dirs['consensus']}/{self.prefix}.masked.dict")
        self.cmds.append(f"{self.gatk} HaplotypeCaller "
                        f"-R {self.dict_dirs['consensus']}/{self.prefix}.masked.fa "
                        f"-I {self.dict_dirs['map']}/{self.prefix}.bam "
                        f"-O {self.dict_dirs['variant']}/{self.prefix}.gatk.vcf")
        self.cmds.append(f"{self.python} {self.dir_bin}/gatk_filter.py "
                        f"{self.dict_dirs['variant']}/{self.prefix}.gatk.vcf "
                        f"{self.dict_dirs['variant']}/{self.prefix}.gatk_filtered.vcf "
                        f"{self.dict_dirs['variant']}/{self.prefix}.gatk_final.vcf")
        self.cmds.append(f"ln -sf {self.dict_dirs['variant']}/{self.prefix}.gatk_final.vcf "
                        f"{self.dict_dirs['variant']}/{self.prefix}.raw.vcf")

    def variant(self):
        logging.info("检测变异")
        self.cmds.append("# 检测变异")
        self.cmds.append(
            f"{self.bedtools} maskfasta -fi {self.bwaidx} "
            f"-bed {self.dict_dirs['map']}/{self.prefix}.lowcovmask.bed "
            f"-fo {self.dict_dirs['consensus']}/{self.prefix}.masked.fa"
        )
        # gatk流程需要.fai和.dict
        self.cmds.append(f"{self.samtools} faidx {self.dict_dirs['consensus']}/{self.prefix}.masked.fa")
        self.freebayes_variants() # 一分钟的力量, 弃用
        # self.gatk_variant_pipe()
        self.cmds.append(f"""
#靶向区域筛选
{self.bedtools} intersect -header -wa \
    -a {self.dict_dirs['variant']}/{self.prefix}.raw.vcf \
    -b {self.dir_etc}/wy_long_amplicon.bed \
    > {self.dict_dirs['variant']}/{self.prefix}.target.vcf
#snpEff注释
{self.snpEff} NC_045512.2 \
    -htmlStats {self.dict_dirs['variant']}/{self.prefix}.snpEff.html \
    -csvStats {self.dict_dirs['variant']}/{self.prefix}.snpEff.csv \
    {self.dict_dirs['variant']}/{self.prefix}.target.vcf \
    > {self.dict_dirs['variant']}/{self.prefix}.snpEff.vcf
ln -sf {self.prefix}.snpEff.vcf {self.dict_dirs['variant']}/{self.prefix}.vcf
{self.bgzip} -c {self.dict_dirs['variant']}/{self.prefix}.vcf > {self.dict_dirs['variant']}/{self.prefix}.vcf.gz
{self.bcftools} index {self.dict_dirs['variant']}/{self.prefix}.vcf.gz
#VCF转TABLE表格
{self.python} {self.dir_bin}/vcf2table.py {self.dict_dirs['variant']}/{self.prefix}.vcf \
    {self.dict_dirs['variant']}/{self.prefix}.trans.tsv
{self.vcftools} --remove-indels --recode --recode-INFO-all \
    --vcf {self.dict_dirs['variant']}/{self.prefix}.target.vcf \
    --out {self.dict_dirs['variant']}/{self.prefix}.snps
{self.bgzip} -c {self.dict_dirs['variant']}/{self.prefix}.snps.recode.vcf \
    > {self.dict_dirs['variant']}/{self.prefix}.snps.recode.vcf.gz
{self.bcftools} index {self.dict_dirs['variant']}/{self.prefix}.snps.recode.vcf.gz
        """)

    def consensus(self):
        logging.info("一致性序列")
        self.cmds.append("# 一致性序列")
        self.cmds.append(
            f"{self.bcftools} consensus -p {self.prefix} --haplotype A "
            f"-f {self.dict_dirs['consensus']}/{self.prefix}.masked.fa "
            f"-o {self.dict_dirs['consensus']}/{self.prefix}.consensus.fa "
            f"{self.dict_dirs['variant']}/{self.prefix}.vcf.gz"
        )
        self.cmds.append(f"sed -i 's/>.*/>{self.prefix}/g' {self.dict_dirs['consensus']}/{self.prefix}.consensus.fa")

    def single_tree(self):
        logging.info("单样本进化树")
        self.cmds.append("# 单样本进化树")
        self.cmds.append(f"cat {self.dict_dirs['consensus']}/{self.prefix}.consensus.fa "
                        f"{self.dir_etc}/SARS-CoV-2_isolates13.fasta "
                        f"> {self.dict_dirs['phylogenetic']}/{self.prefix}.isolates13.fasta")
        self.cmds.append(f"{self.mafft} --thread 2 --auto --maxiterate 1000 "
                        f"{self.dict_dirs['phylogenetic']}/{self.prefix}.isolates13.fasta "
                        f"> {self.dict_dirs['phylogenetic']}/{self.prefix}.aln.fa")
        self.cmds.append(f"{self.fasttree} -nt {self.dict_dirs['phylogenetic']}/{self.prefix}.aln.fa "
                        f"> {self.dict_dirs['phylogenetic']}/{self.prefix}.tre")
        self.cmds.append(f"{self.Rscript} {self.dir_bin}/tree4SC2pipe.R {self.dict_dirs['phylogenetic']}/{self.prefix}.tre")
        self.cmds.append(f"{self.python} {self.dir_bin}/magick.py {self.magick} {self.dict_dirs['phylogenetic']}")

    def global_sars2(self):
        logging.info("新冠全球库")
        if self.gsars2:
            self.cmds.append("# 新冠全球库")
            self.cmds.append(f"{self.python} {self.dir_bin}/sars_cov2_global.py "
                            f"-i {self.dict_dirs['consensus']}/{self.prefix}.consensus.fa "
                            f"-o {self.dict_dirs['global_database']} -p {self.prefix}")
            self.cmds.append(f"rm {self.dict_dirs['global_database']}/pipe.sh")

    def remove_intermedia(self):
        logging.info("删除中间文件")
        self.cmds.append("# 删除中间文件")
        self.cmds.extend([
            f"rm {self.dict_dirs['qc']}/{self.prefix}.clean.[12].fq",
            f"rm {self.dict_dirs['map']}/{self.prefix}.*.bam*"
        ])

    def make_result_dirs(self):
        logging.info("创建结果文件目录")
        # if os.path.isdir(self.outdir): # 先别删
        #     shutil.rmtree(self.outdir)
        self.dict_dirs = {}  # {"qc": "1.qc" ...}
        # !顺序不要变
        dirs = ["1.qc", "2.map", "3.variant", "4.consensus", "5.phylogenetic"] # [220628 update] 一分钟的力量, 弃用 
        if self.gsars2: # 如果要做新冠全球库
            dirs.append("6.global_database") 
        for dr in dirs:
            os.makedirs(f"{self.outdir}/{dr}", exist_ok=True)
            _key = re.sub(r"\d\.", "", dr)
            self.dict_dirs[_key] = f"{self.outdir}/{dr}"

    def execute(self):
        """执行全流程"""
        logging.info("执行全流程")
        self.make_result_dirs()
        self.qc()
        self.alignment()
        self.variant()
        self.consensus()
        self.single_tree()    # [220628 update] 一分钟的力量, 弃用 
        self.global_sars2()
        #self.remove_intermedia()
        # 全流程脚本
        with open(f"{self.outdir}/pipe.sh", "wt", encoding="utf-8", newline="") as gh:
            for cmd in self.cmds:
                gh.write(cmd + "\n")
        if not self.dryrun: # 是否执行
            run(f"bash {self.outdir}/pipe.sh", shell=True)

if __name__ == "__main__":
    args = get_args()
    pipe = Amplicon(args)
    pipe.execute()
