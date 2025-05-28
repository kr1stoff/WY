#!/usr/bin/env python

import os
import sys
import re
import gzip
import json
import logging
from pathlib import Path
from collections import namedtuple
from subprocess import run
from Bio import SeqIO
import yaml


def wraprun(cml):
    """包装run减少重复操作"""
    logging.debug(f"运行命令: {cml}")
    res = run(cml, shell=True, capture_output=True, encoding="utf-8")
    return res


class MakeDB():
    """创建Swissprot数据库,包含blastp和diamond库"""
    def __init__(self, infa, outdir) -> None:
        self.infa = infa
        self.outdir = outdir

    def execute(self):
        self.check_params()
        self.get_params()
        self.global_params()
        self.my_makedirs()
        self.generate_taxon_table()
        self.generate_lineages()
        self.divisions_meta()
        self.diamond_blastp_makedb()

    def my_makedirs(self):
        self.dir_divisions = f"{self.outdir}/divisions"
        os.makedirs(self.dir_divisions, exist_ok=True)

    def check_params(self):
        """检查输入参数"""
        if not Path(self.infa).exists():
            sys.exit("输入FA不存在!")

    def global_params(self):
        """一般的全局变量放在这"""
        self.metakdm = ["Archaea", "Bacteria", "Viruses", "Fungi"]

    def get_params(self):
        """获取软件/数据库等参数"""
        yml = Path(__file__).parents[1].joinpath("conf/swissprot_db.yaml")
        dic = yaml.safe_load(open(yml, "rt"))
        Params = namedtuple("Params", dic.keys())
        self.params = Params(**dic)

    def generate_taxon_table(self):
        """生成 taxonomy lineage 表格, 顺便生成注释文件(JSON格式)"""
        dict_anno = {}
        pattern = re.compile("sp\|(.*?) (.*?) OS=(.*?) OX=(\d+) GN=(.*?) PE=(\d*) SV=(\d*)", re.S)
        # taxon 表格
        with open(f"{self.outdir}/uniprot_sprot.taxon2", "wt", encoding="utf-8", newline="") as gh:
            for rec in SeqIO.parse(gzip.open(self.infa, "rt"), "fasta"):
                res = re.findall(pattern, rec.description)
                if res:
                    # ID|登录名, 蛋白名称, 来源物种, taxid, 基因名, 可靠性, 版本号
                    id_acc, name, organism, taxid, gene, existence, ver = res[0]
                    gh.write(f"{rec.id}\t{taxid}\n") # id - taxid 表
                    dict_anno[id_acc] = {}
                    dict_anno[id_acc]["name"] = name
                    dict_anno[id_acc]["organism"] = organism
                    dict_anno[id_acc]["taxid"] = taxid
                    dict_anno[id_acc]["gene"] = gene
                    dict_anno[id_acc]["existence"] = existence
                    dict_anno[id_acc]["ver"] = ver
        # 注释信息
        with open(f"{self.outdir}/uniprot_sprot.anno.json", "wt", encoding='utf-8') as fh:
            fh.write(json.dumps(dict_anno, indent=4, ensure_ascii=False))

    def generate_lineages(self):
        """使用taxonkit生成所有条目的lineage"""
        # linux 命令行操作 
        cml = f"""
set -e
cut -f2 {self.outdir}/uniprot_sprot.taxon2 \
    | {self.params.taxonkit} lineage \
    > {self.outdir}/uniprot_sprot.lineage
paste {self.outdir}/uniprot_sprot.taxon2 {self.outdir}/uniprot_sprot.lineage \
    > {self.outdir}/uniprot_sprot.taxon3
rm {self.outdir}/uniprot_sprot.taxon2 {self.outdir}/uniprot_sprot.lineage
        """
        wraprun(cml)

    def divisions_meta(self):
        """拆分FASTA数据库到字库,再合并成META库."""
        if Path(f"{self.outdir}/meta.fasta").exists():
            os.remove(f"{self.outdir}/meta.fasta")
        # 宏基因组四类
        for kdm in self.metakdm:
            cml = f"""
set -e
grep "{kdm}" {self.outdir}/uniprot_sprot.taxon3 \
    | cut -f1 > {self.dir_divisions}/{kdm}.ids
{self.params.seqtk} subseq {self.infa} {self.dir_divisions}/{kdm}.ids \
    > {self.dir_divisions}/{kdm}.fasta
cat {self.dir_divisions}/{kdm}.fasta >> {self.outdir}/meta.fasta
            """
            wraprun(cml)

    def diamond_blastp_makedb(self):
        """创建diamond数据库"""
        fastas = [f"{self.dir_divisions}/{kdm}.fasta" for kdm in self.metakdm]
        fastas.append(f"{self.outdir}/meta.fasta")
        for fa in fastas:
            cml = f"""
set -e
{self.params.diamond} makedb --in {fa} -d {fa}
{self.params.makeblastdb} -in {fa} -dbtype prot -parse_seqids -out {fa}
            """
            wraprun(cml)
