#!/usr/bin/env python
# 230112 删除GC, 添加平均长度和Q20

import os
import sys
from pathlib import Path
import json
from collections import defaultdict


def repfunc(dict_filter_summery):
    if 'read2_mean_length' in dict_filter_summery:
        mealen = str(int((dict_filter_summery['read1_mean_length']+dict_filter_summery['read1_mean_length'])/2))
    else:
        mealen   = str(dict_filter_summery['read1_mean_length'])
    return mealen

def dict_fastp_outtable(out_dict, fastp_dict):
    # before filtering
    out_dict["before"]["reads"]     = format(fastp_dict["summary"]["before_filtering"]["total_reads"], ",")
    out_dict["before"]["bases"]     = format(fastp_dict["summary"]["before_filtering"]["total_bases"], ",")
    out_dict["before"]["Q20"]       = "{:.2%}".format(fastp_dict["summary"]["before_filtering"]["q20_rate"])
    out_dict["before"]["Q30"]       = "{:.2%}".format(fastp_dict["summary"]["before_filtering"]["q30_rate"])
    out_dict["before"]["meanlen"]    = repfunc(fastp_dict["summary"]["before_filtering"])
    # after filtering

    out_dict["after"]["reads"]     = format(fastp_dict["summary"]["after_filtering"]["total_reads"], ",")
    out_dict["after"]["bases"]     = format(fastp_dict["summary"]["after_filtering"]["total_bases"], ",")
    out_dict["after"]["Q20"]       = "{:.2%}".format(fastp_dict["summary"]["after_filtering"]["q20_rate"])
    out_dict["after"]["Q30"]       = "{:.2%}".format(fastp_dict["summary"]["after_filtering"]["q30_rate"])
    out_dict["after"]["meanlen"]    = repfunc(fastp_dict["summary"]["after_filtering"])
    return out_dict, fastp_dict

def fastp_json2table(fpjson, fptable):
    """
    Translate fastp .json to qc table.
    Input: fastp.json
    Output: fastq_stats.txt
    """
    fastp_json, out_file = fpjson, fptable
    fastp_dict = json.load(open(fastp_json))
    out_dict = defaultdict(dict)
    #提取信息
    out_dict, fastp_dict = dict_fastp_outtable(out_dict, fastp_dict)
    # out
    header_list = ["状态", "总序列数", "总碱基数", '序列平均长度', "Q20", "Q30"]
    ba_dict = {"before": "过滤前", "after": "过滤后"}
    outheaderlist = ["reads", "bases", 'meanlen', "Q20", "Q30"]
    with open(out_file, "wt", encoding="utf-8", newline="") as g:
        g.write("\t".join(header_list) + "\n")
        for ab in ["before", "after"]:
            outlist = [out_dict[ab][ohl] for ohl in outheaderlist]
            outlist.insert(0, ba_dict[ab])
            g.write("\t".join(outlist) + "\n")

if __name__ == "__main__":
    usage = f"Usage:\n\t{Path(sys.argv[0]).name} <fastp_json> <fastp_table>\n"
    if len(sys.argv) != 3:
        print(usage)
        sys.exit(1)
    fastp_json, fastp_table = sys.argv[1:]
    fastp_json2table(fastp_json, fastp_table)
