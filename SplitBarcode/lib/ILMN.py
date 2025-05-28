#!/usr/bin/env python

import sys
import os
import logging
import yaml
from subprocess import run
from glob import glob
import xml.dom.minidom
from Bio.Seq import Seq # third-party


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler)


class ILMN():
    def __init__(self, inyaml):
        self.inyaml = inyaml
        self.dir_conf = os.path.join(sys.path[0], "conf")
        self.dir_template = os.path.join(sys.path[0], "template")
        self.assign_base_info()
        self.assign_soft_dict()
        self.assign_parameters()
        self.judge_read_type()

    def assign_base_info(self):
        logging.info("读取输入配置文件")
        with open(self.inyaml, "rt") as fh:
            self.dict_input = yaml.safe_load(fh)
        self.dir_bcl = self.dict_input["BCL_folder"]
        self.library = self.dict_input["library"]
        self.dir_fastq = self.dict_input["FASTQ_folder"]
        self.outdir = os.path.join(self.dir_fastq, self.library)
        self.mismatch = self.dict_input["mismatch"]
        self.index_type = self.dict_input["index_type"]
        self.dict_samples = self.dict_input["samples"]

    def assign_soft_dict(self):
        """软件路径"""
        logging.info("读取软件配置文件")
        with open(os.path.join(self.dir_conf, "software.yaml"), "rt") as fh:
            self.dict_soft = yaml.load(fh, Loader=yaml.SafeLoader)
        self.bcl2fastq = self.dict_soft["bcl2fastq"]
        self.splitBarcode = self.dict_soft["splitBarcode"]

    def assign_parameters(self):
        """软件配置相关类参数"""
        logging.info("读取软件参数配置文件")
        with open(os.path.join(self.dir_conf, "parameters.yaml"), "rt") as fh:
            self.dict_params = yaml.load(fh, Loader=yaml.SafeLoader)
        self.thread = self.dict_params["thread"]

    def judge_read_type(self):
        """判断测序类型: PE/SE"""
        if os.path.isfile(f"{self.dir_bcl}/RunParameters.xml"):
            dom = xml.dom.minidom.parse(f"{self.dir_bcl}/RunParameters.xml")
            root = dom.documentElement
            IsPairedEnd = root.getElementsByTagName("IsPairedEnd")
            if IsPairedEnd[0].firstChild.data:
                self.read_type = "PE"
            else:
                self.read_type = "SE"
        else:
            logging.error("不存在 RunParameters.xml 文件, BCL文件夹不完整!")

    def generate_samplesheet(self):
        logging.info("制作 samplesheet")
        template = os.path.join(self.dir_template, "VenusTemplate.csv")
        self.samplesheet = f"{self.outdir}/{self.library}.samplesheet.csv"
        _flag = False
        with open(template, "rt") as fh, open(self.samplesheet, "wt", encoding="utf-8", newline="") as gh:
            for line in fh:
                if line.startswith("Sample_ID"):
                    llist = line.strip().split(",")
                    lenum = enumerate(llist)
                    # Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,
                    # I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
                    dict_data_head = {item[1]:item[0] for item in lenum} # [data] 表头
                gh.write(line)
            for name in self.dict_samples:
                outlist = [""] * 12 # samplesheet [data] 12列
                outlist[dict_data_head["Sample_ID"]] = name
                outlist[dict_data_head["index"]] = self.dict_samples[name]["index"]
                if self.index_type == "Dual" and self.dict_input["index2_direction"] == "forward": # 双index,正向
                    outlist[dict_data_head["index2"]] = self.dict_samples[name]["index2"]
                elif self.index_type == "Dual" and self.dict_input["index2_direction"] == "reverse": # 双index,反向
                    outlist[dict_data_head["index2"]] = str(Seq(self.dict_samples[name]["index2"]).reverse_complement())
                gh.write(",".join(outlist) + "\n")

    def mybcl2fastq(self):
        logging.info("bcl2fastq 数据拆分")
        cmd = f"{self.bcl2fastq} --no-lane-splitting --barcode-mismatches {self.mismatch} --processing-threads {self.thread} "\
            f"--runfolder-dir {self.dir_bcl} --input-dir {self.dir_bcl}/Data/Intensities/BaseCalls "\
            f"--sample-sheet {self.samplesheet} --output-dir {self.outdir} "\
            f"> {self.outdir}/logs/bcl2fastq.out 2> {self.outdir}/logs/bcl2fastq.err"
        logging.debug(f"cmd: {cmd}")
        run(cmd, shell=True)

    def find_fastq(self, pattern_fastq):
        """找FASTQ文件"""
        res = glob(pattern_fastq)
        if len(res) == 1:
            fastq = res[0]
        else:
            fastq = "" # 没找到
        return fastq

    def generate_id_fastq(self):
        logging.info("生成 ID-FASTQ 对应表格")
        if os.path.isdir(f"{self.outdir}/Stats"): # 拆分成功标识
            gh = open(f"{self.outdir}/logs/id_fastq.tsv", "wt", encoding="utf-8", newline="")
            for name in self.dict_samples:
                outlist = [name]
                fq1 = self.find_fastq(f"{self.outdir}/{name}*R1*.fastq.gz")
                outlist.append(fq1)
                if self.read_type == "PE": # PE
                    fq2 = self.find_fastq(f"{self.outdir}/{name}*R2*.fastq.gz")
                    outlist.append(fq2)
                gh.write("\t".join(outlist) + "\n")
            gh.close()
        else:
            logging.error("拆分失败! 检查记录文件.")
        
    def mymakedirs(self):
        logging.info("创建目录")
        os.makedirs(f"{self.outdir}/logs", exist_ok=True)
    
    def execute(self):
        self.mymakedirs()
        self.generate_samplesheet()
        self.mybcl2fastq()
        self.generate_id_fastq()
