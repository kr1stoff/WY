#!/usr/bin/env python
import sys
import os
import re
import logging
from pathlib import Path
from openpyxl import Workbook


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

def printUsage():
    logging.error("PROG <lineage_report>")
    sys.exit(1)

if len(sys.argv) != 2:
    logging.error("需要一个参数!")
    printUsage()
elif not os.path.isfile(sys.argv[1]):
    logging.error("Pangolin结果文件不存在!")
    printUsage()
else:
    lineage_report = sys.argv[1]
    dir_lineage = os.path.dirname(lineage_report)
    lineage_report_trans = f"{dir_lineage}/lineage_report_trans.xlsx"
    lineage_report_trans_tab = f"{dir_lineage}/lineage_report_trans.xls"

# 配置文件
HOME = Path(__file__).parents[1]
who_name_cn_file = os.path.join(HOME, "etc/WHOLabel_CN.txt")
cov_lineage_db = os.path.join(HOME, "etc/Cov_Lineage.txt")

# [20220126 update] WHO中文名
who_name_cn_dict = dict()
with open(who_name_cn_file) as fh:
    for line in fh:
        linelist = line.replace("\n", "").split("\t")
        who_name_cn_dict[linelist[0]] = linelist[1]

# 重置表格
remain_list = ["Lineage",
                "Most common countries",
                "Earliest date",
                "Description"]
# IO
cov_lineage_dict = dict()
with open(cov_lineage_db, "rt") as fh:
    header = next(fh)
    header_list = header.strip().split("\t")
    header_enum = enumerate(header_list)
    header_dict = {enum[1]: enum[0] for enum in header_enum}
    for line in fh:
        linelist = line.strip().split("\t")
        lineage = linelist[header_dict["Lineage"]]
        outlist = [linelist[header_dict[rl]] for rl in remain_list]
        cov_lineage_dict[lineage] = outlist

# [20220214 update lanlei] 生成excel格式结果
wb = Workbook()
ws = wb.active
with open(lineage_report, "rt") as fh, open(lineage_report_trans_tab, "wt", newline="", encoding="utf-8") as gh:
    outheader = ["样本名", "WHO命名", "Pangolin谱系", "最普遍国家", "最早发现日期", "描述"]
    gh.write("\t".join(outheader) + "\n")
    ws.append(outheader)
    header = next(fh)
    header_list = header.strip().split(",")
    header_enum = enumerate(header_list)
    header_dict = {enum[1]: enum[0] for enum in header_enum}
    for line in fh:
        linelist        = line.strip().split(",")
        sample          = linelist[header_dict["taxon"]]
        who_name_tmp    = linelist[header_dict["scorpio_call"]]
        # 格式化 scorpio_call 到 WHO label
        if " " in who_name_tmp: # Omicron
            who_name_en  = who_name_tmp.replace("Probable ", "").strip().split(" ")[0]
            if who_name_en in who_name_cn_dict:
                who_name = "{} ({})".format(who_name_en, who_name_cn_dict[who_name_en])
            else:
                who_name = who_name_en
        elif who_name_tmp == "":    # 没有分型
            who_name = "-"
        else:   # 没中文的分型
            who_name = who_name_tmp 
        # 删除参考行 reference NC_045512.2 | OM648778.1
        pattern = re.compile("NC_045512.2|MN908947.3", re.I)
        res = re.findall(pattern, sample) 
        if len(res) != 0: 
            continue
        lineage = linelist[header_dict["lineage"]]
        try:
            outlist = [sample, who_name] + cov_lineage_dict[lineage]
        except Exception:
            logging.debug(lineage + " not in cov_lineage_dict")
            outlist = [sample, who_name] + [lineage, "-", "-", "-"]
        ws.append(outlist)
        gh.write("\t".join(outlist) + "\n")
        # 单样本pangolin结果
        with open(f"{dir_lineage}/{sample}.lineage.tsv", "wt", encoding="utf-8", newline="") as dh:
            dh.write("\t".join(outheader[1:]) + "\n")
            dh.write("\t".join(outlist[1:]))
wb.save(lineage_report_trans)
