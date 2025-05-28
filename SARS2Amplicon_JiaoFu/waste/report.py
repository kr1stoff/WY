#!/usr/bin/env python

import os
import sys
import glob
import zipfile
from multiprocessing import Pool


def base_quantity_tmp(read:int=1):
    ## find fastqc results
    read_pattern = "{0}/1.qc/fastqc/*[Rr_]{1}*fastqc/Images/per_base_quality.png".\
    format(single_result_dir, read)
    #~ catch exception
    try:
        read_base_quantity_abs = glob.glob(read_pattern)[0]
    except IndexError or Exception:
        raise Exception ("ERROR - report.py - Not find {}".format(read_pattern))
    read_base_quantity_relative = read_base_quantity_abs.replace(single_result_dir + os.sep, "")
    return read_base_quantity_relative

def copy_report_module():
    """
    cp 
    1. report_module, report.pl and src/
    2. cov_linage, snp tree, multiseq tree
    /sdbb/bioinfor/mengxf/Project/1.microbial_WGS/results/20220125-SC2NGS-T9E9S9T9S/PUMCH496
    """
    # [unuseful] 改版之后没有了项目名和样本名替换
    with open(report_module_pl) as f, open(report_module_pl_copy, "wt", encoding="utf-8", newline="") as g:
        text = f.read().lstrip()
        outtext_base_quantity = text.replace("{{:read1_base_quantity:}}", read1_base_quantity).\
            replace("{{:read2_base_quantity:}}", read2_base_quantity)
        outtext = outtext_base_quantity.replace("{{:项目名:}}", taskid).replace("{{:样本名:}}", samp).strip()
        g.write(outtext + "\n")
    ## linux commands
    #~ check src/, if exists, delete it!
    if os.path.isdir(single_src):
        os.system("rm -r {}".format(single_src))
    os.system("cp -r {} {}".format(report_module_src, single_result_dir))
    #~ single sample cov lineage
    os.makedirs(single_src_source)
    cmd = """
        grep -E '样本名|{0}' {1}/cov_lineage/lineage_report_trans.xls | \
            cut -f2- \
            > {2}/src/source/lineage_single.tsv
        """.format(samp, library_result_dir, single_result_dir)
    os.system(cmd)

def generate_report():
    #~ symlink
    snp_tree_link   = single_src_source + os.sep + "snp_tree_rect.png"
    mseq_tree_link  = single_src_source + os.sep + "mseq_tree_rect.png"
    lineage_trans   = single_src_source + os.sep + "lineage_report_trans.xls"
    os.symlink("../../../snp_phylogenetic_tree/rectangular.png", snp_tree_link)
    os.symlink("../../../multiseq_phylogenetic_tree/rectangular.png", mseq_tree_link)
    os.symlink("../../../cov_lineage/lineage_report_trans.xls", lineage_trans)
    #~ generate index.html
    os.chdir(single_result_dir)
    os.system("/usr/bin/perl " + report_module_pl_copy)
    os.remove(report_module_pl_copy)
    os.chdir(current_work_dir)

def zip_dir(dirpath):
    """
    parameter - dirpath:        target directory path.
    execute   - dirpath.zip:    cropress target directory to .zip
    """
    out_fullname = dirpath + ".zip"
    zip = zipfile.ZipFile(out_fullname, "w", zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(dirpath):
        arcpath = root.replace(os.path.dirname(dirpath), "")
        for file in files:
            zip.write(os.path.join(root, file), arcname=os.path.join(arcpath, file))
    zip.close()


# report module position
program_current_dir     = sys.path[0]
current_work_dir        = os.getcwd()
report_module_dir       = program_current_dir + os.sep + "report_module"
report_module_pl        = report_module_dir + os.sep + "report.pl"
report_module_src       = report_module_dir + os.sep + "src"
sars2_env               = os.path.split(program_current_dir)[0] + os.sep + "libs/sars_cov2.env"

# __main__
assert len(sys.argv) == 2, "ERROR - report.py - Just need 1 Arguments, <TASK ID>!"
taskid = sys.argv[1]
# find result directory
with open(sars2_env) as fh:
    for line in fh:
        if "EARTH_ANALYSIS" in line:
            earth_analysis = line.strip().split("=")[1].replace("\"", "")
# library position
library_result_dir  = "{}/{}".format(earth_analysis, taskid)
# main
single_result_dir_list = list()
sample_path_table   = library_result_dir + os.sep + "intermedia/sample_path_table.txt"
with open(sample_path_table) as fh:
    for line in fh:
        samp = line.strip().split("\t")[0]
        single_result_dir           = "{}/{}".format(library_result_dir, samp)
        report_module_pl_copy       = single_result_dir + os.sep + "report.pl"
        single_src                  = single_result_dir + os.sep + "src"
        single_src_source           = single_src + os.sep + "source"
        single_result_dir_list.append(single_result_dir)
        #` find fastqc base quantity graph
        read1_base_quantity = base_quantity_tmp(1)
        read2_base_quantity = base_quantity_tmp(2)
        #` main process
        copy_report_module()
        generate_report()

# multiprocessing zip
with Pool(12) as ph:
    ph.map(zip_dir, single_result_dir_list)
    
# [20220325] logs directories to logs.zip AnalysisResults/1438/logs --> Analysis/1438/1438_log.zip
logs_dir = library_result_dir + os.sep + "logs"
log_zip = "{}/{}_log.zip".format(library_result_dir.replace("Results",""), taskid)
zip_dir(logs_dir)
os.system("mv {}.zip {}".format(logs_dir, log_zip))
