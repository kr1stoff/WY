#!/usr/bin/env python

import os
import sys
import logging
import re
from pathlib import Path
import pdb
# custom
from file.ei1010 import common, myfastq


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)


class EI1010:
    """
    工作流脚本处理类
        indict      任务-样本名字典
        task_id     任务id
        cursor      MySQL启用中的游标
    """

    def __init__(self, indict, cursor):
        self.indict = indict
        self.cursor = cursor
        self.task_id = self.indict["task_id"]
        self.dir_analysis = f"/sdbb/Earth/Analysis/{self.task_id}"
        self.dir_results = f"/sdbb/Earth/AnalysisResults/{self.task_id}"
        self.log = f"/sdbb/Earth/Analysis/{self.task_id}/polling.log"
        self.samples = self.indict["samples"].keys()  # 样本number列表
        self.popen_returncode = 0  # 默认有异常,运行成功后更新
        self.params = common.get_params()
        # 参考序列
        sql_cmd = "SELECT * FROM weiyuan.tb_dict_ref_gene"
        self.cursor.execute(sql_cmd)
        refs = self.cursor.fetchall()
        self.dict_ref = {str(ref["ref_id"]): ref["ref_seq"] for ref in refs}
        # 进化树
        sql_cmd = "SELECT * FROM weiyuan.tb_dict_evolutionary_tree"
        self.cursor.execute(sql_cmd)
        refs = self.cursor.fetchall()
        self.dict_tree = {str(ref["tree_id"]): ref["tree_ref_fa"] for ref in refs}

    def prepare_bed(self):
        """准备bed文件"""
        in_excel = self.indict["primer_seq"]
        out_bed = common.excel2bed(in_excel, f"{self.dir_analysis}/primers.bed")
        return out_bed

    def my_get_infastq(self, outdict):
        """获取样本表fastq, 写入outdict字典"""
        for name in self.indict["samples"]:
            outdict["samples"][name] = []
            outdict["samples"][name].append(self.indict["samples"][name]["sequencing_data1"])
            if self.indict["samples"][name]["sequencing_data2"] != "":  # 双端
                outdict["samples"][name].append(self.indict["samples"][name]["sequencing_data2"])

    def my_get_infastq_idseqcomplete(self, outdict):
        """[IDseqComplete专用] 修剪FASTQ至100bp. 写入outdict字典"""

        def repfunc1(fq):
            # 重复步骤fq1/2修剪
            sffx_pttrn = re.compile('.fq.gz|.fastq.gz|.fq|.fastq')
            sffx_100 = '.fastq.gz' if Path(fq).name.endswith('.gz') else '.fastq'
            # [230309 BUG] "100bp fastq.gz"输入后缀.fastq会报错
            fq_100_bsnm = re.sub(sffx_pttrn, sffx_100, Path(fq).name)
            fq_100 = f'{self.dir_analysis}/rawdata/{fq_100_bsnm}'
            cmd = myfastq.AnyTo100BP(fq, fq_100).any_to_100bp_command()  # 修剪fq, 命令
            outdict["samples"][name].append(fq_100)
            cmds.append(cmd)

        # [230222 zhangkai] FASTQ切成100bp self.my_get_infastq(outdict)
        Path(f'{self.dir_analysis}/rawdata').mkdir(parents=True, exist_ok=True)
        cmds = []
        for name in self.indict["samples"]:
            logging.debug(name)
            fq1 = self.indict["samples"][name]["sequencing_data1"]
            fq2 = self.indict["samples"][name]["sequencing_data2"]
            outdict["samples"][name] = []
            repfunc1(fq1)
            if fq2 != "":  # 双端
                repfunc1(fq2)
        common.commands2shell(cmds, f'{self.dir_analysis}/rawdata/batch_trimfq.sh')
        os.system(f'cat {self.dir_analysis}/rawdata/batch_trimfq.sh | {self.params.parallel} -j 12')

    def my_get_infasta(self, outdict):
        """获取样本表fasta, 写入outdict字典"""
        for name in self.indict["samples"]:
            outdict["samples"][name] = self.indict["samples"][name]["fasta"]

    def my_get_tree_ref(self, outdict):
        """
        获取进化树参考序列, 写入outdict字典. 
        1.进化树参考是文件夹,里面有很多FA; 2.可以是一个文件夹
        """
        tree_ids = self.indict["tree_data"].split(",")
        for tid in tree_ids:
            dir_tree_ref = self.dict_tree[tid]
            for rfile in os.listdir(dir_tree_ref):
                outdict["samples"][os.path.basename(rfile)] = os.path.join(dir_tree_ref, rfile)

    def my_get_ref(self, outdict):
        """获取参考序列,写入outdict"""
        ref_id = self.indict["reference_genome"]
        ref_fa = self.dict_ref[ref_id]
        outdict["reference"] = ref_fa  # 通用流程不做核酸蛋白注释

    def my_get_pathogen_type(self, outdict):
        """获取病原类型,写入outdict"""
        dict_pathogen = {
            "病毒": "virus",
            "细菌": "bacteria",
            "真菌": "fungi"
        }
        outdict["kindom"] = dict_pathogen[self.indict["pathogen_type"]]

    def trace_tree_shell(self, tree_method, outdict):
        """给trace_tree准备shell"""
        inyaml = f"{self.dir_analysis}/trace_tree.yaml"
        common.dict2yaml(outdict, inyaml)
        # referenced before assignment
        # WGS
        cml = f"python3 /sdbb/share/pipeline/TracePathogen/main.py -p wgs -i {inyaml}"
        if tree_method in ["SNP-FA", "SNP-FQ"]:
            cml = f"python3 /sdbb/share/pipeline/TracePathogen/main.py -p snp -i {inyaml}"
        elif tree_method == "CORE":
            cml = f"python3 /sdbb/share/pipeline/TracePathogen/main.py -p core -i {inyaml}"
        return cml

    def amplicon_upload_return(self):
        """扩增子流程上传和返回MySQL是一样的, 重复步骤"""
        if self.popen_returncode == 0:
            common.lunx_upload(self.dir_analysis, self.dir_results, task_id=self.task_id)
            common.return2mysql(self.task_id, self.samples, self.cursor, return_fasta=True)
        else:
            logging.info(f"popen_returncode: {self.popen_returncode}")
        # [240223 ypy] AnalysisResult 结果可删除, 更改权限777
        os.system(f"chmod -R 777 /sdbb/Earth/AnalysisResults/{self.task_id}")

    def ngs_sars2_wga(self):
        """二代测序新冠全基因组扩增子分析路程"""
        # 准备
        dict_sars_wga = {
            "result_dir": "/sdbb/Earth/Analysis",
            "library": str(self.task_id),
            "samples": {},
            "global_sars2": False,
            "bed": ""
        }
        out_bed = self.prepare_bed()
        if out_bed != "":
            dict_sars_wga["bed"] = out_bed  # bed提交了吗
        if self.indict["sars_cov_2_database"] == "是":
            dict_sars_wga["global_sars2"] = True  # 新冠全球库做不做
        self.my_get_infastq(dict_sars_wga)
        yaml_sars_wga = f"{self.dir_analysis}/ngs_sars2_wga.yaml"
        common.dict2yaml(dict_sars_wga, yaml_sars_wga)
        # 运行
        cml = f"python /sdbb/share/pipeline/Sars_Cov2_Amplicon/main.py -i {yaml_sars_wga}"
        logging.info(f"cml: {cml}")
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()

    def ngs_sars2_wga_i·dseqcomplete(self):
        """IDseqComplete新冠试剂盒专用流程"""
        dict_sars_wga = {
            "result_dir": "/sdbb/Earth/Analysis",
            "library": str(self.task_id),
            "samples": {},
            "global_sars2": False,
            "bed": "/sdbb/share/pipeline/Sars_Cov2_Amplicon/etc/SARS-COV-2.FLEX_With_AddonV1V2O_primer_info.tab"
        }
        if self.indict["sars_cov_2_database"] == "是":
            dict_sars_wga["global_sars2"] = True  # 新冠全球库做不做
        # [230222 zhangkai] FASTQ切成100bp self.my_get_infastq(dict_sars_wga)
        self.my_get_infastq_idseqcomplete(dict_sars_wga)
        yaml_sars_wga = f"{self.dir_analysis}/ngs_sars2_wga.yaml"
        common.dict2yaml(dict_sars_wga, yaml_sars_wga)
        # 运行
        cml = f"python /sdbb/share/pipeline/Sars_Cov2_Amplicon/main.py -i {yaml_sars_wga} --trim_software ktrim"
        logging.info(f"cml: {cml}")
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()

    def ngs_sars_wastewater_surveillance(self):
        """[230411] 二代测序新冠污水监测流程"""
        cml = f"python /sdbb/share/pipeline/SARS2Wastewater/main.py "\
            f"-i {self.dir_analysis}/task_samples.yaml -o /sdbb/Earth/Analysis"
        logging.info(cml)
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        os.system(f'mv {self.dir_analysis}/.log {self.dir_analysis}/logs')
        self.amplicon_upload_return()

    def ont_sars_wastewater_surveillance(self):
        """[230425] 三代测序新冠污水监测流程"""
        cml = f"python /sdbb/share/pipeline/SARS2Wastewater/main.py " \
            f"-i {self.dir_analysis}/task_samples.yaml -o /sdbb/Earth/Analysis -p ONT"
        logging.info(cml)
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        os.system(f'mv {self.dir_analysis}/.log {self.dir_analysis}/logs')
        self.amplicon_upload_return()

    def ngs_amplicon_general(self):
        """二代测序通用扩增子分析流程"""
        dict_amplicon_general = {
            "result_dir": "/sdbb/Earth/Analysis",
            "library": str(self.task_id),
            "samples": {},
            "bed": ""
        }
        # 获取ref
        self.my_get_ref(dict_amplicon_general)  # 通用流程不做核酸蛋白注释
        # bed
        out_bed = self.prepare_bed()
        if out_bed != "":
            dict_amplicon_general["bed"] = out_bed  # bed提交了吗
        dict_amplicon_general["trace"] = False  # 通用流程不做进化树
        dict_amplicon_general["kindom"] = "virus"
        self.my_get_infastq(dict_amplicon_general)
        yaml_amplicon_general = f"{self.dir_analysis}/ngs_amplicon_general.yaml"
        common.dict2yaml(dict_amplicon_general, yaml_amplicon_general)
        # 运行
        cml = f"python /sdbb/share/pipeline/AmpliconGP/main.py -i {yaml_amplicon_general}"
        logging.info(f"cml: {cml}")
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()

    def trace_tree(self):
        """溯源进化树, 全基因组、核心基因、SNP"""
        """二代测序通用扩增子分析流程"""
        dict_trace_tree = {
            "result_dir": "/sdbb/Earth/Analysis",
            "library": str(self.task_id),
            "samples": {}
        }
        # WGS全长 1.没有进化树参考 2.有进化树参考
        if self.indict["tree_method"] == "全长序列" and self.indict["tree_data"] == "":
            self.my_get_infasta(dict_trace_tree)
            cml = self.trace_tree_shell("WGS", dict_trace_tree)
        elif self.indict["tree_method"] == "全长序列" and self.indict["tree_data"] != "":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_tree_ref(dict_trace_tree)
            cml = self.trace_tree_shell("WGS", dict_trace_tree)
        # 基因片段，默认和全基因组一样
        elif self.indict["sequence_type"] == "基因片段":
            self.my_get_infasta(dict_trace_tree)
            cml = self.trace_tree_shell("WGS", dict_trace_tree)
        # 核心基因 细菌真菌 1.没有进化树参考 2.有进化树参考
        elif self.indict["tree_method"] == "coregene" and self.indict["tree_data"] == "":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            cml = self.trace_tree_shell("CORE", dict_trace_tree)
        elif self.indict["tree_method"] == "coregene" and self.indict["tree_data"] != "":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            self.my_get_tree_ref(dict_trace_tree)
            cml = self.trace_tree_shell("CORE", dict_trace_tree)
        # SNP 1.fasta 2.fastq
        elif self.indict["tree_method"] == "SNP" and self.indict["input_fq_or_fa"] == "fasta":
            self.my_get_infasta(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            self.my_get_ref(dict_trace_tree)
            cml = self.trace_tree_shell("SNP-FA", dict_trace_tree)
        elif self.indict["tree_method"] == "SNP" and self.indict["input_fq_or_fa"] == "fastq":
            self.my_get_infastq(dict_trace_tree)
            self.my_get_pathogen_type(dict_trace_tree)
            self.my_get_ref(dict_trace_tree)
            cml = self.trace_tree_shell("SNP-FQ", dict_trace_tree)
        else:
            logging.error(f"不存在的进化树参数组合!!! sequence_type:{self.indict['sequence_type']}, "
                          f"tree_method:{self.indict['tree_method']}, tree_data: {self.indict['tree_data']}")
            return None
        # 运行
        logging.info(f"cml: {cml}")
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        if self.popen_returncode == 0:
            common.lunx_upload(self.dir_analysis, self.dir_results, task_id=self.task_id, batch=True)
            common.return2mysql_batch(self.task_id, self.samples, self.cursor)
        else:
            logging.info(f"popen_returncode: {self.popen_returncode}")
        # [240223 ypy] AnalysisResult 结果可删除, 更改权限777
        os.system(f"chmod -R 777 /sdbb/Earth/AnalysisResults/{self.task_id}")

    def ont_virus_assembly(self):
        """[240228] ONT病毒无参组装"""
        cml = f"python /sdbb/share/pipeline/ONTwholeGenome/main.py " \
            f"-i {self.dir_analysis}/task_samples.yaml -o /sdbb/Earth/Analysis -w virassy"
        logging.info(cml)
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()

    def ont_bacteria_assembly(self):
        """[240229] ONT细菌无参组装"""
        cml = f"python /sdbb/share/pipeline/ONTwholeGenome/main.py " \
            f"-i {self.dir_analysis}/task_samples.yaml -o /sdbb/Earth/Analysis -w batassy"
        logging.info(cml)
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()

    def hyb_virus_assembly(self):
        """[240304] HYB病毒无参组装"""
        cml = f"python /sdbb/share/pipeline/HYBwholeGenome/main.py " \
            f"-i {self.dir_analysis}/task_samples.yaml -o /sdbb/Earth/Analysis -w virassy"
        logging.info(cml)
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()

    def hyb_bacteria_assembly(self):
        """[240304] HYB细菌无参组装"""
        cml = f"python /sdbb/share/pipeline/HYBwholeGenome/main.py " \
            f"-i {self.dir_analysis}/task_samples.yaml -o /sdbb/Earth/Analysis"
        logging.info(cml)
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()

    def ont_virus_map(self):
        """[240228] ONT病毒有参组装"""
        cml = f"python /sdbb/share/pipeline/ONTwholeGenome/main.py " \
            f"-i {self.dir_analysis}/task_samples.yaml -o /sdbb/Earth/Analysis -w virmap"
        logging.info(cml)
        self.popen_returncode = common.wrap_popen(cml, self.task_id, self.cursor, self.log)
        self.amplicon_upload_return()
