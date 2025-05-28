#!/usr/bin/env python
import sys
import os
import glob
import re
import json
import pymysql
import time
import shutil
from pprint import pprint
import pdb


######## CMD Options ########
# read library id from command line
def printUsage():
    print("PROG <library_id>\nNeed 1 Argument!")
    sys.exit(1)
## Argument number
if len(sys.argv) != 2:
    printUsage()
## argv1 == integer?
try:
    int(sys.argv[1])
except ValueError:
    printUsage()
except Exception:
    printUsage()

#===============================================================================
scripts_cache               = "/sdbb/bioinfor/mengxf/Project/1.microbial_WGS/scripts/cache"
#===============================================================================

setting_rawfq_id = 1
settint_work_id = 2
library_id = sys.argv[1]

# library dict
library_dict = dict()

# connect database hospital_lims
connection = pymysql.connect(host="127.0.0.1", 
                            user="hospital_lims",
                            password="tZkALtMsEXPrKz7Y",
                            database="hospital_lims",
                            charset='utf8mb4',
                            cursorclass=pymysql.cursors.DictCursor)
# select from database hospital_lims
with connection.cursor() as cursor:
    ## table "setting" select: 1.setting_bcl_dir 2.setting_work_dir
    ### 1
    sql = "SELECT `url` FROM `setting` WHERE `id`= {}".format(setting_rawfq_id)
    cursor.execute(sql)
    result = cursor.fetchone()
    ### {'url': '/home/id_seq/ID-seq/Share'}
    setting_bcl_dir = result["url"]
    ### 2
    sql = "SELECT `url` FROM `setting` WHERE `id`= {}".format(settint_work_id)
    cursor.execute(sql)
    result = cursor.fetchone()
    ### {'url': '/home/id_seq/WORK'}
    setting_work_dir = result["url"]
    
    ## table "library"
    ### name:library_name, hear:chip_id, operator:operator
    sql = "SELECT `name`, `hear`, `operator` FROM `library` WHERE `id`= {}".format(library_id)
    cursor.execute(sql)
    result = cursor.fetchone()
    ### {'hear': '000H3KCTJ', 'name': '20211211-SC2NGS', 'operator': '九号操作人'}
    library_dict["library_name"] = result["name"]
    library_dict["chip_id"]      = result["hear"]
    library_dict["operator"]     = result["operator"]
    sql = "SELECT `fastq1`, `fastq2`, `sublibrary_number` FROM `sublibrary` WHERE `library_id`= {}".format(library_id)
    cursor.execute(sql)
    results = cursor.fetchall()
    library_dict["samples"] = dict()
    # [{'fastq1': 'PUMCH495_S3_R1_001.fastq.gz','fastq2': 'PUMCH495_S3_R2_001.fastq.gz','sublibrary_number': '211231_PUMCH495'},
    #  {'fastq1': 'PUMCH496_S4_R1_001.fastq.gz','fastq2': 'PUMCH496_S4_R2_001.fastq.gz','sublibrary_number': '211231_PUMCH496'}]
    for item in results:
        sample_name = item["sublibrary_number"]
        library_dict["samples"][sample_name]                = dict()
        library_dict["samples"][sample_name]["sample_name"] = sample_name
        library_dict["samples"][sample_name]["fq1_base"]    = item["fastq1"]
        library_dict["samples"][sample_name]["fq2_base"]    = item["fastq2"]
# close connection
connection.close()

# library dict part
search_dir_list = glob.glob("{}/*{}".format(setting_bcl_dir, library_dict["chip_id"]))
## raw fastq dir: 1.exists, 2.unique
if len(search_dir_list) == 0 or not os.path.isdir(search_dir_list[0]):
    raise Exception ("ERROR - read_mysql.py - Runfolder dir don't exists!")
elif len(search_dir_list) != 1:
    raise Exception ("ERROR - read_mysql.py - CHIP ID not a unique match!")
else:
    ### BaseCalls {runfolder_dir}/Data/Intensities/BaseCalls
    library_dict["rawdata_dir"] = search_dir_list[0]
# results dir
# /home/id_seq/WORK/IDseqV2/20211209-SC2NGS-TEST09
library_dict["result_dir"]   = "{}/IDseqV2/{}".format(setting_work_dir, library_dict["library_name"])

# write json
library_json_file = library_dict["rawdata_dir"] + os.sep + "library.json"
with open(library_json_file, "wt", encoding='utf-8', newline="") as g:
    json.dump(library_dict, g, indent=4, ensure_ascii=False)

# write library_id for idseq-AS pipe
library_id_file = library_dict["rawdata_dir"] + os.sep + "library_id"
with open(library_id_file, "wt", newline="", encoding="utf-8") as g:
    g.write(library_id + "\n")

# write work.env
workenv_file = scripts_cache + os.sep + "work.env"
if os.path.isfile(workenv_file):
    os.remove(workenv_file)
with open(workenv_file, "wt", newline="", encoding="utf-8") as g:
    g.write("export LIBRARY_NAME={}\n".format(library_dict["library_name"]))
    g.write("export RESULT_DIR={}\n".format(library_dict["result_dir"]))
    g.write("export RAWDATA_DIR={}\n".format(library_dict["rawdata_dir"]))
    g.write("export HOME=/home/id_seq\n")
current_time = time.strftime("%Y%m%d_%H%M%S", time.localtime())
workenv_file_bak = workenv_file + "_" + current_time
shutil.copy(workenv_file, workenv_file_bak)
