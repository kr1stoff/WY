#!/usr/bin/env python

import os
import sys
import yaml
import logging
import shutil
from subprocess import run
from Bio import SeqIO
from pathlib import Path
# custom
from lib import common


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def cml2shell(cmds, shell):
    """
    命令list写到shell脚本里
        cmds:  命令列表
        shell: 输出脚本
    """
    with open(shell, "wt", encoding="utf-8", newline="") as gh:
        for cmd in cmds:
            gh.write(cmd + "\n")


class SARS2:
    def __init__(self, args):
        self.home = Path(__file__).resolve().parents[1]
        self.mylog = common.MyLoggingInfo()
        self.cmds = list()
        self.trim_software = args.trim_software

        # basic
        with open(args.inyaml, "rt") as fh:
            self.dict_input = yaml.safe_load(fh)
        libr = self.dict_input["library"]
        self.outdir = self.dict_input["result_dir"] + os.sep + libr
        self.dict_sample = self.dict_input["samples"]
        self.global_sars2 = self.dict_input["global_sars2"]
        self.bed = self.dict_input["bed"]

        # software
        with open(self.home.joinpath("conf/software.yaml"), "rt") as fh:
            dict_soft = yaml.load(fh, Loader=yaml.SafeLoader)
        self.python = dict_soft["python"]
        self.bwa = dict_soft["bwa"]
        self.snpEff = dict_soft["snpEff"]
        self.bcftools = dict_soft["bcftools"]
        self.fasttree = dict_soft["fasttree"]
        self.Rscript = dict_soft["Rscript"]
        self.mafft = dict_soft["mafft"]
        self.magick = dict_soft["magick"]
        self.perl = dict_soft["perl"]
        self.pangolin = dict_soft["pangolin"]

        # parameters
        with open(self.home.joinpath("conf/parameters.yaml"), "rt") as fh:
            dict_params = yaml.load(fh, Loader=yaml.SafeLoader)
        self.thread = dict_params["thread"]
        self.parallel_number = dict_params["parallel_number"]
        with open(self.home.joinpath("conf/databases.yaml"), "rt") as fh:
            dict_db = yaml.safe_load(fh)
            self.bwadb_prefix = dict_db["bwadb_prefix"]

    def link_data(self):
        self.mylog.info("软链接原始数据")
        fq_suff = '.fastq'
        for name in self.dict_sample:
            # [240222] 所有FASTQ, 压缩和未压缩需要一致
            fq_suff = ".fastq.gz" if self.dict_sample[name][0].endswith(".gz") else ".fastq"
            common.link_exist(self.dict_sample[name][0], f"{self.outdir}/rawdata/{name}.1{fq_suff}")
            if len(self.dict_sample[name]) > 1:
                common.link_exist(self.dict_sample[name][1], f"{self.outdir}/rawdata/{name}.2{fq_suff}")
        return fq_suff

    def parallel_variant_anaysis(self, fq_suff):
        self.mylog.info("并行跑单样本基础分析流程")
        cmds_para_var = list()  # 并行分析流程列表
        _bed = f"-t {self.bed}" if self.bed != "" else ""  # 提交了bed文件
        _global_sars2 = "--global_sars2" if self.global_sars2 else ""  # True, 使用新冠全球库
        for name in self.dict_sample:
            if len(self.dict_sample[name]) > 1:  # 双端read加read2 fastq
                _fastq2 = f"-I {self.outdir}/rawdata/{name}.2{fq_suff}"
            else:
                _fastq2 = ""
            cmd = f"""
{self.python} {self.home}/sars_cov2_pipe.py \
    -i {self.outdir}/rawdata/{name}.1{fq_suff} \
    -o {self.outdir}/{name} -p {name} \
    --trim_software {self.trim_software} \
    {_bed} {_global_sars2} {_fastq2} \
    > {self.outdir}/logs/{name}.variant.out 2> {self.outdir}/logs/{name}.variant.err
            """
            cmds_para_var.append(cmd)
        cml2shell(cmds_para_var, f"{self.outdir}/Shell/vars.sh")
        self.cmds.append("# 并行跑单样本基础分析流程")
        self.cmds.append(f"""
{self.python} {self.home}/bin/kparallel.py -p {self.parallel_number} {self.outdir}/Shell/vars.sh \
    > {self.outdir}/logs/vars.out 2> {self.outdir}/logs/vars.err
        """)

    def trace_tree(self):
        """溯源进化树"""
        self.mylog.info("进化树流程")
        # SNP进化树
        os.makedirs(f"{self.outdir}/SNPTree", exist_ok=True)
        snp_cmds = list()
        snp_cmds.append(f"ls {self.outdir}/*/3.variant/*.snps.recode.vcf.gz > {self.outdir}/intermedia/vcf_list.txt")
        snp_cmds.append(
            f"{self.bcftools} merge -m snps -f PASS,. --force-samples --output-type v "
            f"--file-list {self.outdir}/intermedia/vcf_list.txt -o {self.outdir}/SNPTree/merged.vcf"
        )
        snp_cmds.append(f"{self.python} {self.home}/bin/vcf2phylip.py -m 1 -i {self.outdir}/SNPTree/merged.vcf "
                        f"--output-folder {self.outdir}/SNPTree")
        snp_cmds.append(f"{self.fasttree} -nt {self.outdir}/SNPTree/merged.min1.phy "
                        f"> {self.outdir}/SNPTree/merged.min1.tre")
        snp_cmds.append(f"{self.Rscript} {self.home}/bin/tree.R {self.outdir}/SNPTree/merged.min1.tre")
        cml2shell(snp_cmds, f"{self.outdir}/Shell/SNPTree.sh")
        self.cmds.append("# SNP进化树")
        self.cmds.append(f"bash {self.outdir}/Shell/SNPTree.sh > {self.outdir}/logs/SNPTree.out "
                         f"2> {self.outdir}/logs/SNPTree.err")
        self.cmds.append(
            f"{self.python} {self.home}/bin/magick.py {self.magick} {self.outdir}/SNPTree "
            f"> {self.outdir}/logs/SNPTree_magick.out 2> {self.outdir}/logs/SNPTree_magick.err"
        )
        # MSA进化树
        os.makedirs(f"{self.outdir}/MSATree", exist_ok=True)
        msa_cmds = list()
        msa_cmds.append(f"cat {self.bwadb_prefix} {self.outdir}/*/4.consensus/*.consensus.fa "
                        f"> {self.outdir}/MSATree/merged.fa")
        msa_cmds.append(f"{self.mafft} --auto --maxiterate 1000 {self.outdir}/MSATree/merged.fa "
                        f"> {self.outdir}/MSATree/merged.aln.fa")
        msa_cmds.append(f"{self.fasttree} -nt {self.outdir}/MSATree/merged.aln.fa > {self.outdir}/MSATree/merged.tre")
        msa_cmds.append(f"{self.Rscript} {self.home}/bin/tree.R {self.outdir}/MSATree/merged.tre")
        cml2shell(msa_cmds, f"{self.outdir}/Shell/MSATree.sh")
        self.cmds.append("# MSA进化树")
        self.cmds.append(f"bash {self.outdir}/Shell/MSATree.sh > {self.outdir}/logs/MSATree.out "
                         f"2> {self.outdir}/logs/MSATree.err")
        self.cmds.append(f"{self.python} {self.home}/bin/magick.py {self.magick} {self.outdir}/MSATree "
                         f"> {self.outdir}/logs/MSATree_magick.out 2> {self.outdir}/logs/MSATree_magick.err")

    def mypangolin(self):
        self.mylog.info("Pangolin分型")
        os.makedirs(f"{self.outdir}/cov_lineage", exist_ok=True)
        self.cmds.append(f"""
# Pangolin分型
export PATH="$PATH:{os.path.dirname(self.pangolin)}"
cat {self.bwadb_prefix} {self.outdir}/*/4.consensus/*.consensus.fa > {self.outdir}/cov_lineage/merged.fa
echo pangolin时间1: $(date)
{self.pangolin} {self.outdir}/cov_lineage/merged.fa -t {self.thread} -o {self.outdir}/cov_lineage \
    > {self.outdir}/logs/pangolin.out 2> {self.outdir}/logs/pangolin.err
echo pangolin时间2: $(date)
{self.python} {self.home}/bin/cov_lineage_list.py {self.outdir}/cov_lineage/lineage_report.csv \
    > {self.outdir}/logs/cov_lineage.out 2> {self.outdir}/logs/cov_lineage.err
        """)

    def make_quality_summary(self):
        self.mylog.info("质控汇总表")
        self.cmds.append(f"""
# 质控汇总表
{self.python} {self.home}/utils/toolbox.py get-quality-summary -d {self.outdir} \
    2> {self.outdir}/logs/quality_summary.log
        """)

    def copy_upload(self):
        # [220826 update] /sdbb/Earth/Analysis 目录权限关闭，软连接都改成复制
        self.mylog.info("上载结果文件")
        cmds_link_res = list()
        if os.path.isdir(f"{self.outdir}/Upload"):  # 如果存在先删掉,ln目录容易产生嵌套
            shutil.rmtree(f"{self.outdir}/Upload")
        for name in self.dict_sample:
            dir_samp = f"{self.outdir}/{name}"  # 文库结果目录/样本名
            upload_samp = f"{self.outdir}/Upload/{name}"  # 文库结果目录/Upload/样本名
            os.makedirs(f"{upload_samp}/source", exist_ok=True)
            _dirs = ["source", "1.qc", "2.map", "3.variant"]
            for dr in _dirs:
                os.makedirs(f"{upload_samp}/{dr}", exist_ok=True)
            cmds_link_res.append(f"""
# 上载样本: {name}
cp -rf {dir_samp}/1.qc/before \
    {dir_samp}/1.qc/after \
    {dir_samp}/1.qc/{name}.html \
    {dir_samp}/1.qc/{name}.basic.stat.txt \
    {dir_samp}/1.qc/{name}.detail.stat.txt \
    -t {upload_samp}/1.qc
cp -f {dir_samp}/1.qc/before/{name}.1_fastqc/Images/per_base_quality.png \
    {upload_samp}/1.qc/{name}.1_before_per_base_quality.png
cp -f {dir_samp}/1.qc/after/{name}.clean.1_fastqc/Images/per_base_quality.png \
    {upload_samp}/1.qc/{name}.1_after_per_base_quality.png 
cp -f {dir_samp}/2.map/{name}.bam* {upload_samp}/2.map
cp -f {dir_samp}/2.map/{name}.genome_coverage_depth*png {upload_samp}/2.map
cp -f {dir_samp}/3.variant/{name}.trans.tsv {upload_samp}/3.variant/{name}.trans.tsv
cp -f {dir_samp}/3.variant/{name}.trans.xlsx {upload_samp}/3.variant/{name}.trans.xlsx
cp -f {dir_samp}/4.consensus/{name}.consensus.fa {upload_samp}/{name}.fa
cp -rf {self.home}/src {upload_samp}/source
cp -f {self.outdir}/cov_lineage/lineage_report_trans.xls {upload_samp}/source/lineage_report_trans.xls
cp -f {self.outdir}/cov_lineage/lineage_report_trans.xlsx {upload_samp}/source/lineage_report_trans.xlsx
cp -f {self.outdir}/cov_lineage/{name}.lineage.tsv {upload_samp}/source/{name}.lineage.tsv
cp -f {self.outdir}/quality_summary.tsv {upload_samp}
            """)
            if len(self.dict_sample[name]) > 1:  # 双端
                cmds_link_res.append(f"""
cp -f {dir_samp}/1.qc/before/{name}.2_fastqc/Images/per_base_quality.png \
    {upload_samp}/1.qc/{name}.2_before_per_base_quality.png
cp -f {dir_samp}/1.qc/after/{name}.clean.2_fastqc/Images/per_base_quality.png \
    {upload_samp}/1.qc/{name}.2_after_per_base_quality.png
                """)
            if self.global_sars2:  # True, 使用新冠全球库
                cmds_link_res.append(f"cp -rf {dir_samp}/6.global_database {upload_samp}")
        cml2shell(cmds_link_res, f"{self.outdir}/Shell/copy_upload.sh")
        self.cmds.append(f"""
# 上载结果文件
echo copy upload 时间1: $(date)
bash {self.outdir}/Shell/copy_upload.sh > {self.outdir}/logs/copy_upload.out 2> {self.outdir}/logs/copy_upload.err
echo copy upload 时间2: $(date)
        """)

    def report(self):
        self.mylog.info("生成报告")
        rep_cmds = list()
        for name in self.dict_sample:
            rep_cmds.append(f"{self.perl} {self.home}/bin/report.pl {self.outdir}/Upload/{name} {name}")
        cml2shell(rep_cmds, f"{self.outdir}/Shell/report.sh")
        self.cmds.append("# 生成报告")
        self.cmds.append(f"bash {self.outdir}/Shell/report.sh > {self.outdir}/logs/report.out "
                         f"2> {self.outdir}/logs/report.err")

    def zip_results(self):
        self.mylog.info("压缩结果")
        zip_cmds = list()
        zip_cmds.append(f"cd {self.outdir}/Upload")  # linux zip 要切目录
        for name in self.dict_sample:
            zip_cmds.append(f"zip -r {name}.zip {name}")
        cml2shell(zip_cmds, f"{self.outdir}/Shell/zip_dir.sh")
        self.cmds.append("# 压缩结果")
        self.cmds.append(f"bash {self.outdir}/Shell/zip_dir.sh > {self.outdir}/logs/zip_dir.out "
                         f"2> {self.outdir}/logs/zip_dir.err")

    def make_result_dirs(self):
        self.mylog.info("创建结果目录")
        dirs = ["rawdata", "database", "logs", "Upload", "Shell", "intermedia"]
        for dr in dirs:
            os.makedirs(f"{self.outdir}/{dr}", exist_ok=True)

    def execute(self, dryrun):
        self.cmds.append("echo '流程开始时间 --> '$(date)")
        self.make_result_dirs()
        fastq_suffix = self.link_data()
        self.parallel_variant_anaysis(fastq_suffix)
        self.mypangolin()
        self.make_quality_summary()
        self.copy_upload()
        self.report()
        self.zip_results()
        self.cmds.append("echo '流程结束时间 --> '$(date)")
        cml2shell(self.cmds, f"{self.outdir}/Shell/all.sh")
        if not dryrun:
            # 运行
            run(f"time bash {self.outdir}/Shell/all.sh > {self.outdir}/logs/all.out 2> {self.outdir}/logs/all.err",
                shell=True)
