#!/usr/bin/env python
import sys
import os
import glob
import pymysql
from pprint import pprint

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
    
######## Libs ########
index5_file = "/sdbb/bioinfor/mengxf/Project/1.microbial_WGS/scripts/libs/Index_IGT_I5.txt"
index7_file = "/sdbb/bioinfor/mengxf/Project/1.microbial_WGS/scripts/libs/Index_IGT_I7.txt"

# id assign
setting_bcl_id = 1
settint_work_id = 2
library_id = sys.argv[1]    ## example: 62
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
    sql = "SELECT `url` FROM `setting` WHERE `id`= {}".format(setting_bcl_id)
    cursor.execute(sql)
    result = cursor.fetchone()
    ### {'url': '/home/id_seq/ID-seq/Share'}
    setting_bcl_dir = result["url"]
    ### 2
    sql = "SELECT `url` FROM `setting` WHERE `id`= {}".format(settint_work_id)
    cursor.execute(sql)
    result = cursor.fetchone()
    # pprint(result)
    ### {'url': '/home/id_seq/WORK'}
    setting_work_dir = result["url"]
    
    ## table "library"
    ### name:library_name, hear:chip_id, operator:operator
    sql = "SELECT `name`, `hear`, `operator` FROM `library` WHERE `id`= {}".format(library_id)
    cursor.execute(sql)
    result = cursor.fetchone()
    # pprint(result)
    ### {'name': '20211210-SC2NGS-TEST10', 'hear': 'CHIP0000', 'operator': '九号操作人'}
    library_dict["library_name"] = result["name"]
    library_dict["chip_id"]      = result["hear"]
    library_dict["operator"]     = result["operator"]

    ## table "sublibrary"
    ### index:DualIndex, ssublibrary_number:ublibrary_name
    sql = "SELECT `index`, `sublibrary_number` FROM `sublibrary` WHERE `library_id`= {}".format(library_id)
    cursor.execute(sql)
    results = cursor.fetchall()
    # pprint(results)
    ### results
    # [{'index': 'IGT-I5-378#|IGT-I7-90#', 'sublibrary_number': 'C20211210_1|IGT-I5-378#|IGT-I7-90#'}, 
    # {'index': 'IGT-I5-379#|IGT-I7-91#', 'sublibrary_number': 'C20211210_2|IGT-I5-379#|IGT-I7-91#'}, 
    # {'index': 'IGT-I5-380#|IGT-I7-92#', 'sublibrary_number': 'C20211210_3|IGT-I5-380#|IGT-I7-92#'}, 
    # {'index': 'IGT-I5-381#|IGT-I7-93#', 'sublibrary_number': 'C20211210_4|IGT-I5-381#|IGT-I7-93#'}]
    library_dict["samples"] = dict()
    for itemd in results:
        dual_index          = itemd["index"]
        sublibrary_number   = itemd["sublibrary_number"]
        try:
            index5, index7  = dual_index.split("|")
            sample_name     = sublibrary_number.split("|")[0]
        except (ValueError, IndexError):
            raise Exception ("Check index or sublibrary_number, \
                e.g. 'index': 'IGT-I5-378#|IGT-I7-90#', 'sublibrary_number': 'C20211210_1|IGT-I5-378#|IGT-I7-90#'")
        except:
            raise Exception ("Check index or sublibrary_number, \
                e.g. 'index': 'IGT-I5-378#|IGT-I7-90#', 'sublibrary_number': 'C20211210_1|IGT-I5-378#|IGT-I7-90#'")
        library_dict["samples"][sample_name]           = dict()
        library_dict["samples"][sample_name]["index5"] = index5
        library_dict["samples"][sample_name]["index7"] = index7
        ### samples_dict
        # {'C20211210_1': {'index5': 'IGT-I5-378#', 'index7': 'IGT-I7-90#'},
        #  'C20211210_2': {'index5': 'IGT-I5-379#', 'index7': 'IGT-I7-91#'},
        #  'C20211210_3': {'index5': 'IGT-I5-380#', 'index7': 'IGT-I7-92#'},
        #  'C20211210_4': {'index5': 'IGT-I5-381#', 'index7': 'IGT-I7-93#'}}
## close connection
connection.close()

## library dict part
search_dir_list = glob.glob("{}/*{}".format(setting_bcl_dir, library_dict["chip_id"]))
### bcl runfold dir: 1.exists, 2.unique
if len(search_dir_list) != 1:
    raise Exception ("CHIP ID not a unique match!")
elif not os.path.isdir(search_dir_list[0]):
    raise Exception ("Runfolder dir don't exists!")
else:
    library_dict["runfolder_dir"] = search_dir_list[0]
    ## BaseCalls {runfolder_dir}/Data/Intensities/BaseCalls
### rawfastq dir & results dir
# /home/id_seq/WORK/IDseq/20211209-SC2NGS-TEST09
library_dict["rawfastq_dir"] = "{}/IDseq/{}".format(setting_work_dir, library_dict["library_name"])
# /home/id_seq/WORK/IDseqV2/20211209-SC2NGS-TEST09
library_dict["result_dir"]   = "{}/IDseqV2/{}".format(setting_work_dir, library_dict["library_name"])

pprint(library_dict)