#!/usr/bin/env python

import os
import sys
import yaml
import logging
import shutil
from subprocess import run
from Bio import SeqIO
from pathlib import Path
# 自编函数
HOME = Path(__file__).resolve().parents[1]   # "../../" 上上级目录
sys.path.append(str(HOME)) 
from lib import common


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class SARS2():
    def __init__(self, args):
        self.inyaml = args.inyaml
        self.trim_software = args.trim_software
        self.dir_conf = os.path.join(HOME, "conf")
        self.dir_etc = os.path.join(HOME, "etc")
        self.dir_bin = os.path.join(HOME, "bin")
        self.dir_src = os.path.join(HOME, "src")
        self.mylog = common.MyLoggingInfo()
        self.assign_base_info()
        self.assign_soft_dict()
        self.assign_parameters()
        self.cmds = list()

    def assign_base_info(self):
        logging.info("读取输入配置文件")
        with open(self.inyaml, "rt") as fh:
            self.dict_input = yaml.safe_load(fh)
        self.library = self.dict_input["library"]
        self.outdir = self.dict_input["result_dir"] + os.sep + self.library
        self.dict_sample = self.dict_input["samples"]
        self.global_sars2 = self.dict_input["global_sars2"]
        self.bed = self.dict_input["bed"]
        self.company = self.dict_input["company"]

    def assign_soft_dict(self):
        """软件路径"""
        logging.info("读取软件配置文件")
        with open(os.path.join(self.dir_conf, "software.yaml"), "rt") as fh:
            self.dict_soft = yaml.load(fh, Loader=yaml.SafeLoader)
        self.python = self.dict_soft["python"]
        self.bwa = self.dict_soft["bwa"]
        self.snpEff = self.dict_soft["snpEff"]
        self.bcftools = self.dict_soft["bcftools"]
        self.fasttree = self.dict_soft["fasttree"]
        self.Rscript = self.dict_soft["Rscript"]
        self.mafft = self.dict_soft["mafft"]
        self.magick = self.dict_soft["magick"]
        self.perl = self.dict_soft["perl"]
        self.pangolin = self.dict_soft["pangolin"]
        self.parallel = self.dict_soft["parallel"]

    def assign_parameters(self):
        """软件配置相关类参数"""
        logging.info("读取软件参数配置文件")
        with open(os.path.join(self.dir_conf, "parameters.yaml"), "rt") as fh:
            self.dict_params = yaml.load(fh, Loader=yaml.SafeLoader)
        self.thread = self.dict_params["thread"]
        self.parallel_number = self.dict_params["parallel_number"]
        """数据库"""
        with open(os.path.join(self.dir_conf, "databases.yaml"), "rt") as fh:
            self.dict_db = yaml.safe_load(fh)
            self.bwadb_prefix = self.dict_db["bwadb_prefix"]

    def link_data(self):
        def suffix(read:str): # 判断后缀 ".gz"
            return ".fastq.gz" if read.endswith(".gz") else ".fastq"
        self.mylog.info("软链接原始数据")
        for name in self.dict_sample:
            self.read1_suffix = suffix(self.dict_sample[name][0])
            common.link_exist(self.dict_sample[name][0], f"{self.outdir}/rawdata/{name}.1{self.read1_suffix}")
            if len(self.dict_sample[name]) > 1:
                self.read2_suffix = suffix(self.dict_sample[name][1])
                common.link_exist(self.dict_sample[name][1], f"{self.outdir}/rawdata/{name}.2{self.read2_suffix}")

    def parallel_variant_anaysis(self):
        self.mylog.info("并行跑单样本基础分析流程")
        self.para_var_anaysis_shell = f"{self.outdir}/Shell/vars.sh"
        cmds_para_var = list() # 并行分析流程列表
        _bed = f"-t {self.bed}" if self.bed != "" else "" # 提交了bed文件
        _global_sars2 = "--global_sars2" if self.global_sars2 else "" # True, 使用新冠全球库
        for name in self.dict_sample:
            if len(self.dict_sample[name]) > 1: # 双端read加read2 fastq
                _fastq2 = f"-I {self.outdir}/rawdata/{name}.2{self.read2_suffix}"
            else:
                _fastq2 = ""
            cmd = f"""
{self.python} {self.dir_bin}/sars_cov2_pipe.py \
    -i {self.outdir}/rawdata/{name}.1{self.read1_suffix} \
    -o {self.outdir}/{name} -p {name} \
    --trim_software {self.trim_software} \
    {_bed} {_global_sars2} {_fastq2} \
    > {self.outdir}/logs/{name}.variant.out 2> {self.outdir}/logs/{name}.variant.err
            """
            cmds_para_var.append(cmd)
        self.cml2shell(cmds_para_var, f"{self.outdir}/Shell/vars.sh")
        self.cmds.append("# 并行跑单样本基础分析流程")
        self.cmds.append(f"""
{self.python} {self.dir_bin}/kparallel.py -p {self.parallel_number} {self.outdir}/Shell/vars.sh \
    > {self.outdir}/logs/vars.out 2> {self.outdir}/logs/vars.err
        """
        )

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
        snp_cmds.append(f"{self.python} {self.dir_bin}/vcf2phylip.py -m 1 -i {self.outdir}/SNPTree/merged.vcf --output-folder {self.outdir}/SNPTree")
        snp_cmds.append(f"{self.fasttree} -nt {self.outdir}/SNPTree/merged.min1.phy > {self.outdir}/SNPTree/merged.min1.tre")
        snp_cmds.append(f"{self.Rscript} {self.dir_bin}/tree.R {self.outdir}/SNPTree/merged.min1.tre")
        self.cml2shell(snp_cmds, f"{self.outdir}/Shell/SNPTree.sh")
        self.cmds.append("# SNP进化树")
        self.cmds.append(f"bash {self.outdir}/Shell/SNPTree.sh > {self.outdir}/logs/SNPTree.out 2> {self.outdir}/logs/SNPTree.err")
        self.cmds.append(
            f"{self.python} {self.dir_bin}/magick.py {self.magick} {self.outdir}/SNPTree "
            f"> {self.outdir}/logs/SNPTree_magick.out 2> {self.outdir}/logs/SNPTree_magick.err"
        )
        # MSA进化树
        os.makedirs(f"{self.outdir}/MSATree", exist_ok=True)
        msa_cmds = list()
        msa_cmds.append(f"""
cat {self.outdir}/*/4.consensus/*.consensus.fa > {self.outdir}/MSATree/merged.fa
{self.mafft} --auto --maxiterate 1000 {self.outdir}/MSATree/merged.fa > {self.outdir}/MSATree/merged.aln.fa
{self.fasttree} -nt {self.outdir}/MSATree/merged.aln.fa > {self.outdir}/MSATree/merged.tre
{self.Rscript} {self.dir_bin}/tree.R {self.outdir}/MSATree/merged.tre
        """)
        self.cml2shell(msa_cmds, f"{self.outdir}/Shell/MSATree.sh")
        self.cmds.append("# MSA进化树")
        self.cmds.append(f"bash {self.outdir}/Shell/MSATree.sh > {self.outdir}/logs/MSATree.out 2> {self.outdir}/logs/MSATree.err")
        self.cmds.append(f"{self.python} {self.dir_bin}/magick.py {self.magick} {self.outdir}/MSATree "
                        f"> {self.outdir}/logs/MSATree_magick.out 2> {self.outdir}/logs/MSATree_magick.err")

    def mypangolin(self):
        self.mylog.info("Pangolin分型")
        os.makedirs(f"{self.outdir}/cov_lineage", exist_ok=True)
        self.cmds.append(f"""
# Pangolin分型
export PATH="$PATH:{os.path.dirname(self.pangolin)}"
cat {self.outdir}/*/4.consensus/*.consensus.fa > {self.outdir}/cov_lineage/merged.fa
echo pangolin时间1: $(date)
{self.pangolin} {self.outdir}/cov_lineage/merged.fa -t {self.thread} -o {self.outdir}/cov_lineage \
    > {self.outdir}/logs/pangolin.out 2> {self.outdir}/logs/pangolin.err
echo pangolin时间2: $(date)
{self.python} {self.dir_bin}/cov_lineage_list.py {self.outdir}/cov_lineage/lineage_report.csv \
    > {self.outdir}/logs/cov_lineage.out 2> {self.outdir}/logs/cov_lineage.err
        """)

    def cml2shell(self, cmds, shell):
        """
        命令list写到shell脚本里
            cmds:  命令列表
            shell: 输出脚本
        """
        with open(shell, "wt", encoding="utf-8", newline="") as gh:
            for cmd in cmds:
                gh.write(cmd + "\n")

    def copy_upload(self):
        # [220826 update] /sdbb/Earth/Analysis 目录权限关闭，软连接都改成复制
        self.mylog.info("上载结果文件")
        cmds_link_res = list()
        if os.path.isdir(f"{self.outdir}/Upload"): # 如果存在先删掉,ln目录容易产生嵌套
            shutil.rmtree(f"{self.outdir}/Upload")
        for name in self.dict_sample:
            dir_samp = f"{self.outdir}/{name}"              # 文库结果目录/样本名
            upload_samp = f"{self.outdir}/Upload/{name}"    # 文库结果目录/Upload/样本名
            os.makedirs(f"{upload_samp}/source", exist_ok=True)
            _dirs = ["source", "1.qc", "2.map", "3.variant", '5.phylogenetic']
            for dr in _dirs:
                os.makedirs(f"{upload_samp}/{dr}", exist_ok=True)
            cmds_link_res.append(f"""
# 上载样本: {name}
cp -rf {dir_samp}/1.qc/before \
    {dir_samp}/1.qc/after \
    {dir_samp}/1.qc/{name}.html \
    {dir_samp}/1.qc/{name}.basic.stat.txt \
    -t {upload_samp}/1.qc
cp -f {dir_samp}/1.qc/before/{name}.1_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.1_before_per_base_quality.png
cp -f {dir_samp}/1.qc/after/{name}.clean.1_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.1_after_per_base_quality.png
cp -f {dir_samp}/2.map/{name}.bam* {upload_samp}/2.map
cp -f {dir_samp}/2.map/{name}.genome_coverage_depth*png {upload_samp}/2.map
cp -f {dir_samp}/3.variant/{name}.trans.tsv {upload_samp}/3.variant/{name}.trans.tsv
cp -f {dir_samp}/3.variant/{name}.trans.xlsx {upload_samp}/3.variant/{name}.trans.xlsx
cp -f {dir_samp}/4.consensus/{name}.consensus.fa {upload_samp}/{name}.fa
cp -f {dir_samp}/5.phylogenetic/rectangular.png {upload_samp}/5.phylogenetic/rectangular.png
cp -rf {self.outdir}/SNPTree {upload_samp}/source
cp -rf {self.outdir}/MSATree {upload_samp}/source
cp -rf {self.dir_src} {upload_samp}/source
cp -f {self.outdir}/cov_lineage/lineage_report_trans.xls {upload_samp}/source/lineage_report_trans.xls
cp -f {self.outdir}/cov_lineage/lineage_report_trans.xlsx {upload_samp}/source/lineage_report_trans.xlsx
cp -f {self.outdir}/cov_lineage/{name}.lineage.tsv {upload_samp}/source/{name}.lineage.tsv
            """)
            if len(self.dict_sample[name]) > 1: # 双端
                cmds_link_res.append(f"""
cp -f {dir_samp}/1.qc/before/{name}.2_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.2_before_per_base_quality.png
cp -f {dir_samp}/1.qc/after/{name}.clean.2_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.2_after_per_base_quality.png
                """)
            if self.global_sars2:   # True, 使用新冠全球库
                cmds_link_res.append(f"cp -rf {dir_samp}/6.global_database {upload_samp}")
        self.cml2shell(cmds_link_res, f"{self.outdir}/Shell/copy_upload.sh")
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
            # rep_cmds.append(f"{self.perl} {self.dir_bin}/report.pl {self.outdir}/Upload/{name} {name}")
            rep_cmds.append(f"""
{self.Rscript} {self.dir_bin}/report_word.R \
    {self.outdir}/Upload/{name} \
    {self.outdir}/materials/{name}.sample_info.txt \
    {self.company}
            """)
        self.cml2shell(rep_cmds, f"{self.outdir}/Shell/report.sh")
        self.cmds.append(f"""
# 生成报告
cat {self.outdir}/Shell/report.sh | {self.parallel} -j {self.parallel_number} > {self.outdir}/logs/report.out 2> {self.outdir}/logs/report.err
        """)

    def zip_results(self):
        self.mylog.info("压缩结果")
        zip_cmds = list()
        zip_cmds.append(f"cd {self.outdir}/Upload") # linux zip 要切目录
        for name in self.dict_sample:
            zip_cmds.append(f"zip -r {name}.zip {name}")
        self.cml2shell(zip_cmds, f"{self.outdir}/Shell/zip_dir.sh")
        self.cmds.append("# 压缩结果")
        self.cmds.append(f"echo zip_dir 时间1: $(date)") # 一分钟
        self.cmds.append(f"bash {self.outdir}/Shell/zip_dir.sh > {self.outdir}/logs/zip_dir.out 2> {self.outdir}/logs/zip_dir.err")
        self.cmds.append(f"echo zip_dir 时间2: $(date)") # 一分钟

    def make_result_dirs(self):
        self.mylog.info("创建结果目录")
        dirs = ["rawdata", "database", "logs", "Upload", "Shell", "intermedia"]
        for dr in dirs:
            os.makedirs(f"{self.outdir}/{dr}", exist_ok=True)

    def execute(self, dryrun):
        self.cmds.append("echo '流程开始时间 --> '$(date)")
        self.make_result_dirs()
        self.link_data()
        self.parallel_variant_anaysis()
        self.trace_tree()     # [220628 update] 一分钟的力量，弃用 
        self.mypangolin()
        self.copy_upload()
        self.report()
        self.zip_results()
        self.cmds.append("echo '流程结束时间 --> '$(date)")
        self.cml2shell(self.cmds, f"{self.outdir}/Shell/all.sh")
        if not dryrun:
            # 运行
            run(f"time bash {self.outdir}/Shell/all.sh > {self.outdir}/logs/all.out 2> {self.outdir}/logs/all.err", shell=True)
