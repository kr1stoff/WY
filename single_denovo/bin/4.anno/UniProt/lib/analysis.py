#!/usr/bin/env python

import sys
import logging
from pathlib import Path
import json
from openpyxl import Workbook
import pandas as pd
sys.path.append(Path(__file__).resolve().parent)
from . import make_db


class Analysis():
    """
    对输入FA进行 Swissprot 数据库分析
    1.比对软件可选 blastp/diamond
    2.数据库可选 宏/细菌/真菌/病毒
    """
    def __init__(self, infa, softw_align, db_select, db_path, outable):
        self.infa = infa
        self.softw_align = softw_align
        self.db_select = db_select
        self.db_path = db_path
        self.outable = outable
        self.outdir = Path(outable).resolve().parent
        make_db.MakeDB.get_params(self)
  
    def blastp(self, bpath):
        cml = f"""
{self.params.blastp} -num_threads {self.params.threads} -evalue 1e-5 \
    -query {self.infa} -db {bpath} \
    -max_hsps 1 -max_target_seqs 1 \
    -outfmt '6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp' \
    -out {self.outdir}/swissprot_aln.out
        """
        make_db.wraprun(cml)

    def diamond(self, dpath):
        cml = f"""
{self.params.diamond} blastp -p {self.params.threads} -e 1e-5 \
    -q {self.infa} -d {dpath} --max-hsps 1 --max-target-seqs 1 \
    --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
    --out {self.outdir}/swissprot_aln.out
        """
        make_db.wraprun(cml)

    def alignment(self):
        if self.db_select == "Meta":
            db_blastp_path = f"{self.db_path}/meta.fasta"
            db_diamond_path = f"{self.db_path}/meta.fasta.dmnd"
        else:
            db_blastp_path = f"{self.db_path}/divisions/{self.db_select}.fasta"
            db_diamond_path = f"{self.db_path}/divisions/{self.db_select}.fasta.dmnd"
        if self.softw_align.lower() == "diamond":
            self.diamond(db_diamond_path)
        elif self.softw_align.lower() == "blastp":
            self.blastp(db_blastp_path)
        else:
            sys.exit(f"不存在的比对软件: {self.softw_align}")

    def get_annotations(self):
        json_anno = f"{self.db_path}/uniprot_sprot.anno.json"
        self.dict_anno = json.load(open(json_anno))

    def make_result(self):
        """生成最终结果表格tsv和excel"""
        sp_excel = f"{self.outdir}/{Path(self.outable).stem}.xlsx"
        headers_aln = ['qseqid', 'qlen', 'sseqid', 'slen', 'pident', 'length', 'mismatch', 'gapopen', 
                        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovhsp']
        hi = {item[1]:item[0] for item in enumerate(headers_aln)} # 表头下标
        with open(f"{self.outdir}/swissprot_aln.out") as fh:
            headers = ['预测基因ID', 'UniProtID', '蛋白名称', '微生物体', '覆盖度(%)', '置信度(%)', 
                        '预测基因长度', '参考长度', '比对长度', '比对得分']
            df = pd.DataFrame(columns=headers)
            for line in fh:
                llist = line.strip().split("\t")
                uniprot_id = llist[hi["sseqid"]].replace("sp|", "")
                subdict_anno = self.dict_anno[uniprot_id]
                protein_name = subdict_anno["name"]
                organism = subdict_anno["organism"]
                outlist = [llist[hi["qseqid"]], uniprot_id, protein_name, organism, float(llist[hi["qcovhsp"]]), 
                        llist[hi["pident"]], llist[hi["qlen"]], llist[hi["slen"]], int(llist[hi["length"]]), llist[hi["bitscore"]]]
                df = pd.concat([df, pd.DataFrame([outlist], columns=headers)])

        df_out = df[df["覆盖度(%)"] > 80].sort_values(by="比对长度", ascending=False) # 覆盖度>80%, 使用比对长度重排序
        df_out.to_csv(self.outable, sep="\t", index=False, encoding="utf-8")
        df_out.to_excel(sp_excel, index=False, encoding="utf-8")

    def execute(self):
        """执行"""
        self.alignment()
        self.get_annotations()
        self.make_result()
