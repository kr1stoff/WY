#!/usr/bin/env python

import os
import re
import string
import logging
import yaml


# 设置运行日志
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)


class CheckWebInput():
    """
    检查VENUS网页端传入MySQL参数是否符合现有标准
        task_id     项目id
        cursor      MySQL已启动的游标
    """
    def __init__(self, task_id, cursor):
        self.task_id = task_id
        self.cursor = cursor
        self.error_message = ""
        self.dir_conf = os.path.join(os.path.dirname(__file__),"conf")
        self.mysql2indict() # 任务名及关联样本信息大字典
        self.get_pathogen_dict() # 获取病原库信息字典
        self.confirm_workflow()
        self.check_all()

    def update_error_message(self, new_error):
        if new_error not in self.error_message:
            self.error_message += new_error

    def mysql2indict(self):
        """mysql转到任务样本总字典"""
        # 任务表
        sql = f"select * from tb_analyse_task where task_id={self.task_id}"
        self.cursor.execute(sql)
        self.indict = self.cursor.fetchone() # cursor 类型改为了字典 zhuzai
        # 样本表
        sql = f"select * from tb_task_sample where task_id={self.task_id}"
        self.cursor.execute(sql)
        tuple_sampledicts = self.cursor.fetchall() # 字典组成的元组
        self.indict["samples"] = dict()
        for samp in tuple_sampledicts:
            self.indict["samples"][samp["sample_number"]] = samp
        return self.indict

    def get_pathogen_dict(self):
        """获取病原库信息字典"""
        sql_cmd = "SELECT * FROM weiyuan.tb_dict_pathogen"
        self.cursor.execute(sql_cmd)
        dict_pathogen = self.cursor.fetchall()
        self.dict_pathogen = {str(dp["pathogen_list_id"]): dp["pathogen_list_name"] for dp in dict_pathogen} 

    def confirm_workflow(self):
        """新版找流程名,要把字段组合先配置好"""
        # 判断字段: level2_pipeline, pe_or_se, pathogen_list_name, pathogen_typ, reference_genome
        level2_pipeline = self.indict["level2_pipeline"]
        # 病原数据库id转name
        # pdb.set_trace()
        pathogen_id = str(self.indict["pathogen_id"])
        pathogen_name = self.dict_pathogen[pathogen_id] if pathogen_id != "None" else ""
        pe_or_se = self.indict["pe_or_se"]
        pathogen_type = self.indict["pathogen_type"]
        reference_genome = self.indict["reference_genome"]
        fields = "_".join([level2_pipeline, pe_or_se, pathogen_name, pathogen_type, reference_genome])
        # 流程-字段对照表
        with open(os.path.join(self.dir_conf,"pipe_field_table.yaml")) as fh:
            dict_pipe_field = yaml.safe_load(fh)
        for pfk in dict_pipe_field:
            pattern = re.compile(pfk)
            if re.match(pattern, fields):
                logging.debug(f"{pattern}: {fields}")
                self.workflow = dict_pipe_field[pfk]
        if "workflow" not in dir(self): # 迭代后workflow没有被定义
            self.workflow = f"不存在的工作流 <{level2_pipeline} {pe_or_se} {pathogen_name} {pathogen_type}>\n"

    def check_sample_number(self, sample_number):
        """样本编号可以包括:大小写字母,数字以及下划线横线'_-'"""
        sample_number_chars = string.ascii_letters + string.digits + "_-"
        for single in sample_number:
            if single in sample_number_chars:# 不支持中文(u'\u4e00' <= single <= u'\u9fff')
                continue
            else:
                e_message = f"样本编号检查: 限定外特殊字符 <{sample_number}>\n"
                self.update_error_message(e_message)
                break
    
    def is_file(self, infile):
        """输入文件存在且不为目录或其他类型"""
        if not os.path.isfile(infile):
            e_message = f"文件检查: 输入文件不存在或其他类型 <{infile}>\n"
            self.update_error_message(e_message)
    
    def myexists(self, infile):
        if not os.path.exists(infile):
            e_message = f"文件检查: 输入文件不存在或其他类型 <{infile}>\n"
            self.update_error_message(e_message)

    def check_label_compare(self, label_compare, func_flag=0):
        """1.组比较限定字符：大写字母，英文逗号","; 2.存在; 3.符合格式AB,BC..."""
        label_compare_chars = string.ascii_uppercase + ","
        for single in label_compare:
            if single not in label_compare_chars:
                func_flag = 1
        if not re.search(r"[A-Z]{2,}", label_compare):
            func_flag = 1
        if (label_compare == "") or (func_flag == 1):
            e_message = f"组名比较检查: 限定外特殊字符或不符合填写格式 <{label_compare}>\n"
            self.update_error_message(e_message)

    def check_all(self):
        """检查所有"""
        if "不存在" in self.workflow: # 工作流存在吗?
            e_message = self.workflow
            self.update_error_message(e_message)
        else:
            for samp in self.indict["samples"]:
                sample_number = self.indict["samples"][samp]["sample_number"]
                sequencing_data1 = self.indict["samples"][samp]["sequencing_data1"]
                sequencing_data2 = self.indict["samples"][samp]["sequencing_data2"]
                label_compare = self.indict["samples"][samp]["label_compare"]
                # 检查所有样本名
                self.check_sample_number(sample_number)
                # 检查使用FASTQ输入的流程
                if self.indict["pe_or_se"]  == "单端":
                    self.is_file(sequencing_data1)
                elif self.indict["pe_or_se"] == "双端":
                    self.is_file(sequencing_data1)
                    self.is_file(sequencing_data2)
                elif self.indict["pe_or_se"] == 'Nanopore' or self.workflow == "nanopore鉴定":
                    self.myexists(sequencing_data1)
                # 检查16S
                if self.workflow == "16S测序":
                    self.check_label_compare(label_compare=label_compare)
