#!/usr/bin/env python

#!/usr/bin/env python

import sys
import os
import logging
import yaml
from subprocess import run
from glob import glob
import re
from lib.ILMN import ILMN # mine


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler)


class MGI(ILMN):
    def __init__(self, inyaml):
        self.inyaml = inyaml
        self.dir_conf = os.path.join(sys.path[0], "conf")
        self.assign_base_info()
        self.assign_soft_dict()
        self.assign_parameters()
        self.judge_read_type()

    def judge_read_type(self):
        """判断测序类型: PE/SE, 不知道MGI有什么参数文件,只能从输入看"""
        read_len = self.dict_input["reads_length"]
        self.read_type = read_len[:2] # SE50 --> SE

    def generate_barcode_file(self):
        logging.info("制作 barcode file")
        self.barcodefile = f"{self.outdir}/{self.library}.barcode.csv"
        gh = open(self.barcodefile, "wt", encoding="utf-8", newline="")
        for name in self.dict_samples:
            outlist = [name, self.dict_samples[name]["index"]]
            if self.index_type == "Dual": # 双index
                outlist[1] += self.dict_samples[name]["index2"]
            gh.write("\t".join(outlist) + "\n")
        gh.close()

    def get_colrows(self):
        """获取Col/Row参数值"""
        cols, rows = [], []
        self.dir_calfiles = f"{self.dir_bcl}/L01/calFile"
        for cfile in glob(f"{self.dir_calfiles}/*cal"):
            cal = os.path.basename(cfile)
            pattern = re.compile("C(\d+)R(\d+).cal")
            try:
                col, row = re.findall(pattern, cal)[0]
                cols.append(int(col))
                rows.append(int(row))
            except Exception:
                logging.debug(f"calfile: {cal}")
        if cols == [] or rows == []:
            logging.critical(f"没有calFile, 检查目录! <{self.dir_bcl}/L01/calFile>")
            sys.exit(1)
        self.col, self.row = max(cols), max(rows)

    def get_cycle_number(self):
        """获取cycle number参数值"""
        self.cycles = len(glob(f"{self.dir_bcl}/L01/metrics/*QC.csv")) # MGIFBS:-C看metrics下面多少个文件

    def get_barcode_info(self):
        """获取barcode信息"""
        names = list(self.dict_samples.keys())
        bc1_len = len(self.dict_samples[names[0]]["index"]) # barcode1长度
        if self.index_type == "Dual":
            bc2_len = len(self.dict_samples[names[0]]["index2"]) # barcode2长度
            bc2_start = self.cycles - bc2_len + 1  # barcode2开始位置, MGIFBS:测序开始一个bp的碱基矫正,barcode位置+1
            bc1_start = self.cycles - bc2_len - bc1_len + 1 # barcode1开始位置
            self.barcode2_info = f" -i {bc2_start} {bc2_len} {self.mismatch} "
        else:
            bc1_start = self.cycles - bc1_len + 1
            self.barcode2_info = ""
        self.barcode1_info = f" -i {bc1_start} {bc1_len} {self.mismatch} "

    def mysplitbarcode(self):
        logging.info("splitBarcode 数据拆分")
        self.get_colrows()
        self.get_cycle_number()
        self.get_barcode_info()
        cmd = f"{self.splitBarcode} -t {self.thread} -F {self.dir_bcl}/L01/calFile "\
            f"--Col {self.col} --Row {self.row} -B {self.barcodefile} -C {self.cycles} -N Output -o {self.outdir} "
        cmd += (self.barcode1_info + self.barcode2_info)  # -i 52 10 1 -i 62 10 1 # MGIFBS:双index两个-i参数
        cmd += f" > {self.outdir}/logs/splitbarcode.out 2> {self.outdir}/logs/splitbarcode.err"
        logging.debug(f"cmd: {cmd}")
        run(cmd, shell=True)

    def generate_id_fastq(self):
        logging.info("生成 ID-FASTQ 对应表格")
        if os.path.isfile(f"{self.outdir}/Output/L01/SequenceStat.txt"): # 拆分成功标识
            gh = open(f"{self.outdir}/logs/id_fastq.tsv", "wt", encoding="utf-8", newline="")
            for name in self.dict_samples:
                outlist = [name]
                if self.read_type == "PE": # PE
                    fq1 = self.find_fastq(f"{self.outdir}/Output/L01/*{name}*_1.fq.gz")
                    fq2 = self.find_fastq(f"{self.outdir}/Output/L01/*{name}*_2.fq.gz")
                    outlist.extend([fq1, fq2])
                else: # SE,PE read1不一样, 和illumina不一样
                    fq1 = self.find_fastq(f"{self.outdir}/Output/L01/*{name}.fq.gz") # Output/L01/Output_L01_UDB-177.fq.gz
                    outlist.append(fq1)
                gh.write("\t".join(outlist) + "\n")
            gh.close()
        else:
            logging.error(f"拆分失败! 检查记录文件. <{self.outdir}/Output/L01/SequenceStat.txt>")

    def execute(self):
        self.mymakedirs()
        self.generate_barcode_file()
        self.mysplitbarcode()
        self.generate_id_fastq()
