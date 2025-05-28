#!/usr/bin/env python

import os
import sys
import yaml
import logging
from subprocess import run
from Bio import SeqIO
from lib import common

# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler)


class Analysis():
    def __init__(self, inyaml):
        self.inyaml = inyaml
        self.dir_conf = os.path.join(sys.path[0], "conf")
        self.dir_etc = os.path.join(sys.path[0], "etc")
        self.dir_bin = os.path.join(sys.path[0], "bin")
        self.dir_src = os.path.join(sys.path[0], "src")
        self.mylog = common.MyLoggingInfo()
        self.assign_base_info()
        self.assign_soft_dict()
        self.check_input_dbtype()
        self.assign_parameters()
        self.cmds = list()

    def assign_base_info(self):
        logging.info("读取输入配置文件")
        with open(self.inyaml, "rt") as fh:
            self.dict_input = yaml.safe_load(fh)
        self.library = self.dict_input["library"]
        self.outdir = self.dict_input["result_dir"] + os.sep + self.library
        self.dict_sample = self.dict_input['samples']
        self.ref = self.dict_input["reference"]
        self.bed = self.dict_input["bed"]
        self.trace = self.dict_input["trace"]
        if self.trace:
            self.kindom = self.dict_input["kindom"]

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

    def assign_parameters(self):
        """软件配置相关类参数"""
        logging.info("读取软件参数配置文件")
        with open(os.path.join(self.dir_conf, "parameters.yaml"), "rt") as fh:
            self.dict_params = yaml.load(fh, Loader=yaml.SafeLoader)
        self.thread = self.dict_params["thread"]
        self.parallel = self.dict_params["parallel"]

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

    def optimize_gff(self):
        """BCBio_gff的结果不太好用,还要再改动一下"""
        with open(f"{self.outdir}/database/ref/genes.raw.gff", "rt") as fh, \
            open(self.gff, "wt", newline="", encoding="utf-8") as gh:
            for line in fh:
                llist = line.split("\t")
                if line.startswith("#"):
                    gh.write(line)
                elif llist[8] == "\n": # 第9列空值, 解决snpEff gff报错问题
                    llist[8] = f"type={llist[2]}"
                    gh.write("\t".join(llist) + "\n")
                elif "db_xref=" in llist[8] and "gene" in llist[8]: # 修改一下第9列, snpEff用db_xref属性
                    attributes = llist[8].strip().split(";")
                    dict_attr = dict()
                    for attr in attributes:
                        _tmp = attr.split("=")
                        dict_attr[_tmp[0]] = _tmp[1]
                    dict_attr["db_xref"] = dict_attr["gene"]
                    llist[8] = ";".join([f"{k}={v}" for k,v in dict_attr.items()])
                    gh.write("\t".join(llist) + "\n")
                else:
                    gh.write(line)

    def build_snpEff_db(self, records):
        """构建snpEff数据库"""
        from BCBio import GFF
        os.makedirs(f"{self.outdir}/database/ref", exist_ok=True)
        os.makedirs(f"{self.outdir}/database/genomes", exist_ok=True)
        common.link_exist(self.fa, f"{self.outdir}/database/genomes/ref.fa")
        self.gff = f"{self.outdir}/database/ref/genes.gff"
        with open(f"{self.outdir}/database/ref/genes.raw.gff", "wt", encoding="utf-8", newline="") as gh:
            GFF.write(records, gh)
        self.optimize_gff()
        # snpEff.config
        with open(f"{self.dir_etc}/snpeff.config", "rt") as fh, \
            open(f"{self.outdir}/database/snpeff.config", "wt", encoding="utf-8", newline="") as gh:
            for line in fh:
                gh.write(line)
            gh.write("ref.genome : AGP Reference\n")
            fids = [record.id for record in records] # 每个chromosome的id
            gh.write(f"\tref.chromosome : {','.join(fids)}\n")
        # snpEff build
        self.cmds.append("# snpEff建库")
        self.cmds.append(f"snpEff build -c {self.outdir}/database/snpeff.config -dataDir . -gff3 ref "
        f"> {self.outdir}/logs/snpeff_build.out 2> {self.outdir}/logs/snpeff_build.err")

    def check_input_dbtype(self):
        """检查输入的数据库是fasta还是genbank"""
        with open(self.ref, "rt") as fh:
            context = fh.read()
        if "LOCUS" in context:
            self.db_type = "genbank"
        elif ">" in context:
            self.db_type = "fasta"
        else:
            logging.error(f"参考数据库不是FASTA也不是GenBank类型. <{self.ref}>")
            sys.exit(1)

    def prepare_db(self):
        """
        数据库准备, 如果
            1. 输入是genbank格式, 转成fasta格式, bwa建索引
            2. 输入是fasta格式,直接bwa建索引
            3. 其他情况报错
        """
        self.mylog.info("数据库准备")
        self.fa = f"{self.outdir}/database/ref.fa"
        if self.db_type == "genbank": # 变异注释
            self.ref_gbk = self.ref
            common.link_exist(self.ref_gbk, f"{self.outdir}/database/ref.gbk")
            records = [record for record in SeqIO.parse(self.ref_gbk, "genbank")]
            SeqIO.write(records, self.fa, "fasta")
            # snpEff 建库
            self.build_snpEff_db(records)
        elif self.db_type == "fasta": # 变异不注释
            self.ref_fa = self.ref
            # [231229] ypy下载参考是Windows换行符CRLF, 要转一下
            cml = f"""
cp -f {self.ref_fa} {self.fa}
dos2unix {self.fa}
            """
            run(cml, shell=True)
        # 建索引
        self.cmds.append("# 参考数据库建BWA索引")
        self.cmds.append(f"{self.bwa} index {self.outdir}/database/ref.fa "
        f"> {self.outdir}/logs/bwa_index.out 2> {self.outdir}/logs/bwa_index.err")

    def parallel_variant_anaysis(self):
        """并行跑单样本基础分析流程"""
        self.mylog.info("并行跑单样本基础分析流程")
        self.para_var_anaysis_shell = f"{self.outdir}/Shell/vars.sh"
        cmds_para_var = list() # 并行分析流程列表
        tails = " " # 加在基础流程后面, 输入genbank？有没有bed？
        if self.bed != "": tails += f"-t {self.bed} " # 提交了bed文件
        if self.db_type == "genbank": tails += f"--snpEff_config {self.outdir}/database/snpeff.config " # 参考为genbank格式
        for name in self.dict_sample:
            cmd = f"{self.python} {self.dir_bin}/amplicon_sequence_analysis.py "\
            f"-i {self.outdir}/rawdata/{name}.1{self.read1_suffix} -o {self.outdir}/{name} -p {name} "\
            f"--bwa_index {self.outdir}/database/ref.fa "
            cmd += tails # 尾巴加上
            if len(self.dict_sample[name]) > 1: # 双端read加read2 fastq
                cmd += f"-I {self.outdir}/rawdata/{name}.2{self.read2_suffix} "
            cmd += f" > {self.outdir}/logs/{name}.variant.out 2> {self.outdir}/logs/{name}.variant.err " # 加记录文件
            cmds_para_var.append(cmd)
        self.cml2shell(cmds_para_var, f"{self.outdir}/Shell/vars.sh")
        self.cmds.append("# 并行跑单样本基础分析流程")
        self.cmds.append(
            f"{self.python} {self.dir_bin}/kparallel.py -p {self.parallel} {self.outdir}/Shell/vars.sh "
            f"> {self.outdir}/logs/vars.out 2> {self.outdir}/logs/vars.err"
        )
    
    def trace_tree(self):
        """
        溯源进化树(可选)
        病毒:       SNP/WGMSA
        细菌/真菌:  SNP
        """
        self.mylog.info("进化树流程(可选)")
        if self.trace:
            # 如果画进化树,各kindom都做SNP进化树
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
            if self.kindom == "virus": # 全基因组多序列比对进化树WGMSA
                os.makedirs(f"{self.outdir}/MSATree", exist_ok=True)
                msa_cmds = list()
                msa_cmds.append(f"cat {self.fa} {self.outdir}/*/4.consensus/*.consensus.fa > {self.outdir}/MSATree/merged.fa")
                msa_cmds.append(f"{self.mafft} --auto --maxiterate 1000 {self.outdir}/MSATree/merged.fa > {self.outdir}/MSATree/merged.aln.fa")
                msa_cmds.append(f"{self.fasttree} -nt {self.outdir}/MSATree/merged.aln.fa > {self.outdir}/MSATree/merged.tre")
                msa_cmds.append(f"{self.Rscript} {self.dir_bin}/tree.R {self.outdir}/MSATree/merged.tre")
                self.cml2shell(msa_cmds, f"{self.outdir}/Shell/MSATree.sh")
                self.cmds.append("# MSA进化树")
                self.cmds.append(f"bash {self.outdir}/Shell/MSATree.sh > {self.outdir}/logs/MSATree.out 2>{self.outdir}/logs/MSATree.err")
                self.cmds.append(
                    f"{self.python} {self.dir_bin}/magick.py {self.magick} {self.outdir}/MSATree "
                    f"> {self.outdir}/logs/MSATree_magick.out 2> {self.outdir}/logs/MSATree_magick.err"
                )
            elif self.kindom in ["bacteria", "fungi"]:
                pass    # 占个位置
            else:
                logging.error(f"kindom参数错误: kindom:<{self.kindom}>")
        else:
            logging.info(f"不做进化树! trace:<{self.trace}>")

    def cml2shell(self, cmds, shell):
        """
        命令list写到shell脚本里
            cmds:  命令列表
            shell: 输出脚本
        """
        with open(shell, "wt", encoding="utf-8", newline="") as gh:
            for cmd in cmds:
                gh.write(cmd + "\n")

    def link_upload(self):
        self.mylog.info("上载结果文件")
        cmds_link_res = list()
        for name in self.dict_sample:
            dir_samp = f"{self.outdir}/{name}"              # 文库结果目录/样本名
            upload_samp = f"{self.outdir}/Upload/{name}"    # 文库结果目录/Upload/样本名
            os.makedirs(f"{upload_samp}/source", exist_ok=True)
            _dirs = ["source", "1.qc", "2.map", "3.variant", "4.consensus"]
            for dr in _dirs:
                os.makedirs(f"{upload_samp}/{dr}", exist_ok=True)
            cmds_link_res.append(f"""
# 上载样本: {name}
cp -f {dir_samp}/1.qc/before/{name}.1_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.1_before_per_base_quality.png
cp -f {dir_samp}/1.qc/after/{name}.clean.1_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.1_after_per_base_quality.png
cp -rf {dir_samp}/1.qc/before \
    {dir_samp}/1.qc/after \
    {dir_samp}/1.qc/{name}.html \
    {dir_samp}/1.qc/{name}.basic.stat.txt \
    {dir_samp}/1.qc/{name}.detail.stat.txt \
    -t {upload_samp}/1.qc
cp -f {dir_samp}/2.map/{name}.bam \
    {dir_samp}/2.map/{name}.bam_stats.txt \
    {dir_samp}/2.map/{name}.genome_coverage_depth.png \
    {dir_samp}/2.map/{name}.genome_coverage_depth_ylim1000.png \
    -t {upload_samp}/2.map
cp -f {dir_samp}/3.variant/{name}.trans.tsv {upload_samp}/3.variant/{name}.trans.tsv
cp -f {dir_samp}/3.variant/{name}.trans.xlsx {upload_samp}/3.variant/{name}.trans.xlsx
cp -f {dir_samp}/4.consensus/{name}.consensus.fa {upload_samp}/4.consensus/{name}.consensus.fa
cp -rf {self.dir_src} {upload_samp}/source
            """)
            if len(self.dict_sample[name]) > 1: # 双端
                cmds_link_res.append(f"""
cp -f {dir_samp}/1.qc/before/{name}.2_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.2_before_per_base_quality.png
cp -f {dir_samp}/1.qc/after/{name}.clean.2_fastqc/Images/per_base_quality.png {upload_samp}/1.qc/{name}.2_after_per_base_quality.png
                """)
            if self.trace: # 做进化树
                cmds_link_res.append(f"cp -rf {self.outdir}/SNPTree {upload_samp}/source")
                if self.kindom == "virus":
                    cmds_link_res.append(f"cp -rf {self.outdir}/MSATree {upload_samp}/source")
        self.cml2shell(cmds_link_res, f"{self.outdir}/Shell/link_upload.sh")
        self.cmds.append("# 上载结果文件")
        self.cmds.append(f"bash {self.outdir}/Shell/link_upload.sh > {self.outdir}/logs/link_upload.out 2> {self.outdir}/logs/link_upload.err")

    def report(self):
        self.mylog.info("生成报告")
        rep_cmds = list()
        for name in self.dict_sample:
            rep_cmds.append(f"{self.perl} {self.dir_bin}/report.pl {self.outdir}/Upload/{name} {name}")
        self.cml2shell(rep_cmds, f"{self.outdir}/Shell/report.sh")
        self.cmds.append("# 生成报告")
        self.cmds.append(f"bash {self.outdir}/Shell/report.sh > {self.outdir}/logs/report.out 2> {self.outdir}/logs/report.err")

    def zip_results(self):
        self.mylog.info("压缩结果")
        zip_cmds = list()
        zip_cmds.append(f"cd {self.outdir}/Upload") # linux zip 要切目录
        for name in self.dict_sample:
            zip_cmds.append(f"zip -r {name}.zip {name}")
        self.cml2shell(zip_cmds, f"{self.outdir}/Shell/zip_dir.sh")
        self.cmds.append("# 压缩结果")
        self.cmds.append(f"bash {self.outdir}/Shell/zip_dir.sh > {self.outdir}/logs/zip_dir.out 2> {self.outdir}/logs/zip_dir.err")

    def make_result_dirs(self):
        self.mylog.info("创建结果目录")
        dirs = ["rawdata", "database", "logs", "Upload", "Shell", "intermedia"]
        for dr in dirs:
            os.makedirs(f"{self.outdir}/{dr}", exist_ok=True)

    def execute(self, dryrun):
        self.make_result_dirs()
        self.link_data()
        self.prepare_db()
        self.parallel_variant_anaysis()
        self.trace_tree()
        self.link_upload()
        self.report()
        self.zip_results()
        self.cml2shell(self.cmds, f"{self.outdir}/Shell/all.sh")
        if not dryrun:
            # 运行
            run(f"bash {self.outdir}/Shell/all.sh > {self.outdir}/logs/all.out 2> {self.outdir}/logs/all.err", shell=True)
