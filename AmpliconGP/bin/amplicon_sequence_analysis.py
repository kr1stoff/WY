#!/usr/bin/env python

import os
import sys
import shutil
import argparse
import logging
from subprocess import run
import yaml


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)


def get_args():
    parser = argparse.ArgumentParser()
    parser.description="二代测序通用扩增子分析流程, 支持单/双端FASTQ, 可选上传引物BED去引物, 可选上传GBK格式进行snpEff变异注释."
    parser.add_argument("-i", "--fastq1", required=True, help="单端测序输入的FASTQ, 或双端测序的read1.")
    parser.add_argument("-I", "--fastq2", help="双端测序需要填, 双端测序的read.")
    parser.add_argument("-o", "--outdir", default="./analysis", help="结果输出的文件夹.")
    parser.add_argument("-p", "--prefix", default="sample", help="样本名称, 或者输出文件的前缀.")
    parser.add_argument("--bwa_index", dest="bwaidx", required=True, help="BWA索引前缀,并且参数本身也是FASTA.")
    parser.add_argument("--snpEff_config", help="已经建好库的snpEff.config文件路径, 流程自建库使用, 包含ref名目录")
    parser.add_argument("-t", "--bed", help="扩增子的靶向区域, BED文件包含参考名和位置信息.")
    parser.add_argument("-d", "--dryrun", action="store_true", help="是否运行程序, 如选择生成shell脚本不运行.")
    return parser.parse_args()


class Amplicon():
    """
    单样本流程
        args :: 命令行参数 Namespace
    """
    def __init__(self, args):
        self.dir_conf = os.path.join(os.path.dirname(sys.path[0]), "conf")
        self.dir_bin = sys.path[0]
        self.assign_base_info(args)
        self.assign_software()
        self.assign_parameters()
        self.cmds = list()

    def assign_base_info(self, args):
        """读取输入信息"""
        logging.info("读取输入信息")
        self.fastq1 = args.fastq1
        self.outdir = args.outdir
        self.prefix = args.prefix
        self.bwaidx = args.bwaidx
        self.fastq2 = args.fastq2
        self.bed    = args.bed
        self.snpeff_config = args.snpEff_config
        self.dryrun = args.dryrun

    def assign_software(self):
        """软件路径"""
        logging.info("读取软件配置文件")
        with open(os.path.join(self.dir_conf, "software.yaml"), "rt") as fh:
            self.dict_soft = yaml.load(fh, Loader=yaml.SafeLoader)
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

    def assign_parameters(self):
        """软件配置相关类参数"""
        logging.info("读取软件参数配置文件")
        with open(os.path.join(self.dir_conf, "parameters.yaml"), "rt") as fh:
            self.dict_params = yaml.load(fh, Loader=yaml.SafeLoader)
        self.thread = self.dict_params["thread"]
        self.parallel = self.dict_params["parallel"]
        self.min_cov = self.dict_params["min_cov"]
        self.min_alt_count = self.dict_params["min_alt_count"]
        self.min_alt_frac = self.dict_params["min_alt_frac"]
        self.options_fastp = self.dict_params["options_fastp"]

    def qc(self):
        logging.info("原始数据QC")
        if self.fastq2:
            fastq2 = self.fastq2
            fastp2I = f"-I {self.fastq2}"
            fastp2O = f"-O {self.outdir}/1.qc/{self.prefix}.clean.2.fq"
            parse_fastp_mode = "PE"
        else:
            fastq2,fastp2I,fastp2O = "", "", ""
            parse_fastp_mode = "SE"
        self.cmds.append(f"""
# 原始数据QC
mkdir {self.outdir}/1.qc/before {self.outdir}/1.qc/after
{self.fastqc} -t {self.thread} --extract -o {self.outdir}/1.qc/before {self.fastq1} {fastq2} &
{self.fastp} --thread {self.thread} {self.options_fastp} \
    -j {self.outdir}/1.qc/{self.prefix}.json -h {self.outdir}/1.qc/{self.prefix}.html \
    -i {self.fastq1} {fastp2I} \
    -o {self.outdir}/1.qc/{self.prefix}.clean.1.fq {fastp2O}
{self.python} {self.dir_bin}/parse_fastp.json.py \
    -m {parse_fastp_mode} -id {self.prefix} \
    -i {self.outdir}/1.qc/{self.prefix}.json
{self.fastqc} -t {self.thread} --extract \
    {self.outdir}/1.qc/{self.prefix}.clean.1.fq {fastp2O.replace('-O ','')} \
    -o {self.outdir}/1.qc/after
        """)

    def alignment(self):
        logging.info("比对")
        self.cmds.append("# 比对")
        fastq2 = f"{self.outdir}/1.qc/{self.prefix}.clean.2.fq" if self.fastq2 else " "
        self.cmds.append(
            f"{self.bwa} mem -t {self.thread} -M -Y -R '@RG\\tID:{self.prefix}\\tSM:{self.prefix}' "
            f"{self.bwaidx} {self.outdir}/1.qc/{self.prefix}.clean.1.fq {fastq2} | {self.samtools} view -@ 4 -hbS - "
            f"| {self.samtools} sort -@ 4 -o {self.outdir}/2.map/{self.prefix}.sorted.bam -"
        )
        self.cmds.append(f"{self.samtools} index {self.outdir}/2.map/{self.prefix}.sorted.bam")
        # 是否去引物
        if self.bed:
            self.cmds.append("#~ 去引物")
            self.cmds.append(
                f"{self.ivar} trim -q 15 -m 30 -s 4 -e -b {self.bed} "
                f"-p {self.outdir}/2.map/{self.prefix}.ivar_trim "
                f"-i {self.outdir}/2.map/{self.prefix}.sorted.bam"
            )
            self.cmds.append(f"{self.samtools} sort -@ 4 -o {self.outdir}/2.map/{self.prefix}.ivar_trim_sorted.bam {self.outdir}/2.map/{self.prefix}.sorted.bam")
            self.cmds.append(f"cp {self.outdir}/2.map/{self.prefix}.ivar_trim_sorted.bam {self.outdir}/2.map/{self.prefix}.bam")
        else:
            self.cmds.append(f"cp {self.outdir}/2.map/{self.prefix}.sorted.bam {self.outdir}/2.map/{self.prefix}.bam")
        self.cmds.append(f"{self.samtools} index {self.outdir}/2.map/{self.prefix}.bam")
        # 深度统计
        self.cmds.append("#~ 深度统计")
        self.cmds.append(f"{self.samtools} depth -a {self.outdir}/2.map/{self.prefix}.bam > {self.outdir}/2.map/{self.prefix}.bam.depth")
        self.cmds.append(f"{self.Rscript} {self.dir_bin}/genome_coverage.R {self.outdir}/2.map/{self.prefix}.bam.depth {self.outdir}/2.map/{self.prefix}")
        self.cmds.append("#~ 低深度区域")
        self.cmds.append(
            f"{self.bedtools} genomecov -bga -ibam {self.outdir}/2.map/{self.prefix}.bam "
            f"| awk '$4<{self.min_cov}' | {self.bedtools} merge -i - "
            f"| awk '$3-$2>20' 1> {self.outdir}/2.map/{self.prefix}.lowcovmask.bed"
        )

    def variant(self):
        logging.info("检测变异")
        self.cmds.append("# 检测变异")
        self.cmds.append(
            f"{self.bedtools} maskfasta -fi {self.bwaidx} "
            f"-bed {self.outdir}/2.map/{self.prefix}.lowcovmask.bed "
            f"-fo {self.outdir}/4.consensus/{self.prefix}.masked.fa"
        )
        self.cmds.append(
            f"{self.freebayes} -F {self.min_alt_frac} -C {self.min_alt_count} "
            f"-f {self.outdir}/4.consensus/{self.prefix}.masked.fa {self.outdir}/2.map/{self.prefix}.bam "
            f"> {self.outdir}/3.variant/{self.prefix}.raw.vcf"
        )
        # snpeff.config 参数
        if self.snpeff_config: # 有注释
            self.cmds.append(
                f"{self.snpEff} ann -c {self.snpeff_config} -dataDir . "
                f"-htmlStats {self.outdir}/3.variant/{self.prefix}.snpEff.html "
                f"-csvStats {self.outdir}/3.variant/{self.prefix}.snpEff.csv "
                f"ref {self.outdir}/3.variant/{self.prefix}.raw.vcf "
                f"> {self.outdir}/3.variant/{self.prefix}.snpEff.vcf"
            )
            self.cmds.append(f"ln -sf {self.prefix}.snpEff.vcf {self.outdir}/3.variant/{self.prefix}.vcf")
        else: # 没有VCF注释
            self.cmds.append(f"ln -sf {self.prefix}.raw.vcf {self.outdir}/3.variant/{self.prefix}.vcf")
        self.cmds.append(f"{self.bgzip} -c {self.outdir}/3.variant/{self.prefix}.vcf > {self.outdir}/3.variant/{self.prefix}.vcf.gz")
        self.cmds.append(f"{self.bcftools} index {self.outdir}/3.variant/{self.prefix}.vcf.gz")
        # 过滤VCF,并转为table格式
        self.cmds.append(f"{self.python} {self.dir_bin}/vcf2table.py {self.outdir}/3.variant/{self.prefix}.vcf {self.outdir}/3.variant/{self.prefix}.trans.tsv")
        # 挑出snp vcf, 后面建进化树用
        self.cmds.append(
            f"{self.vcftools} --remove-indels --recode --recode-INFO-all "
            f"--vcf {self.outdir}/3.variant/{self.prefix}.raw.vcf "
            f"--out {self.outdir}/3.variant/{self.prefix}.snps"
        )
        self.cmds.append(f"{self.bgzip} -c {self.outdir}/3.variant/{self.prefix}.snps.recode.vcf > {self.outdir}/3.variant/{self.prefix}.snps.recode.vcf.gz")
        self.cmds.append(f"{self.bcftools} index {self.outdir}/3.variant/{self.prefix}.snps.recode.vcf.gz")

    def consensus(self):
        logging.info("一致性序列")
        self.cmds.append("# 一致性序列")
        self.cmds.append(
            f"{self.bcftools} consensus -p {self.prefix} --haplotype A "
            f"-f {self.outdir}/4.consensus/{self.prefix}.masked.fa "
            f"-o {self.outdir}/4.consensus/{self.prefix}.consensus.fa "
            f"{self.outdir}/3.variant/{self.prefix}.vcf.gz"
        )
        self.cmds.append(f"sed -i 's/>.*/>{self.prefix}/g' {self.outdir}/4.consensus/{self.prefix}.consensus.fa")

    def remove_intermedia(self):
        logging.info("删除中间文件")
        self.cmds.append("# 删除中间文件")
        self.cmds.extend([
            f"rm {self.outdir}/1.qc/{self.prefix}.clean.[12].fq",
            f"rm {self.outdir}/2.map/{self.prefix}.*.bam"
        ])

    def make_result_dirs(self):
        logging.info("创建结果文件目录")
        if os.path.isdir(self.outdir):
            shutil.rmtree(self.outdir)
        dirs = ["1.qc", "2.map", "3.variant", "4.consensus"]
        for dr in dirs:
            os.makedirs(f"{self.outdir}/{dr}", exist_ok=True)

    def execute(self):
        """执行全流程"""
        logging.info("执行全流程")
        self.make_result_dirs()
        self.qc()
        self.alignment()
        self.variant()
        self.consensus()
        self.remove_intermedia()
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
