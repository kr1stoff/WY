#!/usr/bin/env python

import os
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
    parser.description="新冠全球库分析脚本。分析路线 'blast --> top n fasta --> MSA'"
    parser.add_argument("-i", "--fasta", required=True, help="输入新冠FASTA序列文件路径")
    parser.add_argument("-o", "--outdir", default="./sc2_global", help="结果输出的文件夹. [default: './sc2_global']")
    parser.add_argument("-p", "--prefix", default="sample", help="样本名称, 或者输出文件的前缀. [default: 'sample']")
    parser.add_argument("-n", "--topN", dest="topn", default="50", 
                        help="比对序列保留的最大序列数,也是进化树选取的代表序列最大数量. [default: 50]")
    parser.add_argument("-d", "--dryrun", action="store_true", help="是否运行程序, 如选择生成shell脚本不运行.")
    return parser.parse_args()


class SarsGlob():
    """新冠全球库分析流程"""
    def __init__(self, args):
        logging.info("读取输入信息")
        self.dir_conf = os.path.abspath(os.path.join(os.path.dirname(__file__), "../conf"))
        self.dir_bin = os.path.dirname(__file__)
        self.assign_base_info(args)
        self.assign_software()
        self.assign_parameters()
        self.cmds = list()

    def assign_base_info(self, args):
        """读取输入信息"""
        logging.info("读取输入信息")
        self.fasta = args.fasta
        self.outdir = args.outdir
        self.prefix = args.prefix
        self.topn = args.topn
        self.dryrun = args.dryrun

    def assign_software(self):
        """软件路径"""
        logging.info("读取软件配置文件")
        with open(os.path.join(self.dir_conf, "software.yaml"), "rt") as fh:
            self.dict_soft = yaml.safe_load(fh)
        self.blastn = self.dict_soft["blastn"]
        self.seqkit = self.dict_soft["seqkit"]
        self.mafft = self.dict_soft["mafft"]
        self.fasttree = self.dict_soft["fasttree"]
        self.Rscript = self.dict_soft["Rscript"]
        self.python = self.dict_soft["python"]
        self.magick = self.dict_soft["magick"]

    def assign_parameters(self):
        """软件配置相关类参数"""
        logging.info("读取软件参数配置文件 & 数据库配置文件")
        with open(os.path.join(self.dir_conf, "parameters.yaml"), "rt") as fh:
            self.dict_params = yaml.safe_load(fh)
        with open(os.path.join(self.dir_conf, "databases.yaml"), "rt") as fh:
            self.dict_db = yaml.safe_load(fh)
        self.thread = self.dict_params["thread"]
        self.db_sars_glob = self.dict_db["sars_global"]

    def blast_top(self):
        logging.info("提取最相似序列")
        self.cmds.append("# 提取最相似序列")
        self.cmds.append(f"{self.blastn} -num_threads {self.thread} -evalue 1e-5  -max_target_seqs {self.topn} "
                        f"-db {self.db_sars_glob} -query {self.fasta} -out {self.outdir}/{self.prefix}.aln "
                        f"-outfmt '7 qseqid sseqid length qlen slen sstart send qstart qend pident nident evalue bitscore'")
        self.cmds.append(f"grep -v '#' {self.outdir}/{self.prefix}.aln | cut -f2 | sort -u "
                        f"> {self.outdir}/{self.prefix}.top{self.topn}.txt")
        self.cmds.append(f"{self.seqkit} grep -j {self.thread} -f {self.outdir}/{self.prefix}.top{self.topn}.txt "
                        f"{self.db_sars_glob} > {self.outdir}/{self.prefix}.top{self.topn}.fa")
        self.cmds.append(f"cat {self.fasta} {self.outdir}/{self.prefix}.top{self.topn}.fa > {self.outdir}/{self.prefix}.top{self.topn}add1.fa")

    def trace_tree(self):
        logging.info("进化树绘制")
        self.cmds.append("# 进化树绘制")
        self.cmds.append(f"{self.mafft} --auto --maxiterate 1000 {self.outdir}/{self.prefix}.top{self.topn}add1.fa > {self.outdir}/{self.prefix}.aln.fa")
        self.cmds.append(f"{self.fasttree} -nt {self.outdir}/{self.prefix}.aln.fa > {self.outdir}/{self.prefix}.tre")
        self.cmds.append(f"{self.Rscript} {self.dir_bin}/tree4global.R {self.outdir}/{self.prefix}.tre")
        self.cmds.append(f"{self.python} {self.dir_bin}/magick.py {self.magick} {self.outdir}")

    def execute(self):
        logging.info("执行全流程")
        os.makedirs(self.outdir, exist_ok=True)
        self.blast_top()
        self.trace_tree()
        # 全流程脚本
        with open(f"{self.outdir}/pipe.sh", "wt", encoding="utf-8", newline="") as gh:
            for cmd in self.cmds:
                gh.write(cmd + "\n")
        if not self.dryrun: # 是否执行
            run(f"bash {self.outdir}/pipe.sh", shell=True)


if __name__ == "__main__":
    args = get_args()
    pipe = SarsGlob(args)
    pipe.execute()
