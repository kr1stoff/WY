#!/usr/bin/env python

##################################################################
## Example:
# python3 ${SCRIPTS_PATH}/utils/pre_process.py library.json
##################################################################

import sys
import os
import json
import glob
import pdb


def printUsage():
    print("PROG <library.json>")
    sys.exit(1)

if len(sys.argv) != 2:
    print("Just need 1 Argument!")
    printUsage()
elif not os.path.isfile(sys.argv[1]):
    print("JSON not exists!")
    printUsage()
else:
    in_json = sys.argv[1]

#~ Open json
with open(in_json) as fh:
    library_dict = json.load(fh)

#~ RawData & Result Path
glb_pattern = "{}/*{}".format(library_dict["rawdata_dir"], library_dict["chip_id"])
glb = glob.glob(glb_pattern)
if len(glb) == 0:
    raise Exception ("ERROR - {} - No Such Directory!".format(glb_pattern))
elif len(glb) > 1:
    raise Exception ("ERROR - {} - No Such Directory!".format(glb_pattern))
else:
    rawdata_dir = glb[0]
result_path = library_dict["result_dir"] + os.sep + library_dict["library_name"]
# pdb.set_trace()
#` intermedia directory
intermedia_path = result_path + os.sep + "intermedia"
if not os.path.isdir(intermedia_path):
    os.makedirs(intermedia_path)
#` sample_path_table.txt
sample_dict = library_dict["samples"]
sample_path_table_file  = intermedia_path + os.sep + "sample_path_table.txt"
outlist = list()
for samp in library_dict["samples"]:
    name  = sample_dict[samp]["sample_name"]
    read1 = rawdata_dir + os.sep + sample_dict[samp]["fq1_base"]
    read2 = rawdata_dir + os.sep + sample_dict[samp]["fq2_base"]
    # check
    if not os.path.isfile(read1):
        raise Exception ("ERROR - {} - Raw FASTQ not exists!".format(read1))
    if not os.path.isfile(read2):
        raise Exception ("ERROR - {} - Raw FASTQ not exists!".format(read2))
    sample_path_line = "{}\t{},{}\n".format(name, read1, read2)
    outlist.append(sample_path_line)
with open(sample_path_table_file, "wt", encoding="utf-8", newline="") as gh:
    for outline in outlist:
        gh.write(outline)
# lookback fasta
lookback_fa = intermedia_path + os.sep + "lookback.fa"
gh = open(lookback_fa, "wt", encoding="utf-8", newline="")
if library_dict["lookback"] == {}:
    os.system("touch " + lookback_fa)
else:
    # "library_name1":  ["A1", "A2", "A3", "..."],
    # "library_name2":  ["A1", "A2", "A3", "..."]
    for lib_name in library_dict["lookback"]:
        for samp in library_dict["lookback"][lib_name]:
            #  AnalysisResults/20211231-SC2NGS-T9E9S9T9S/details/PUMCH496/4.consensus/PUMCH496.consensus.fa
            consensus_fa = "{resultdir}/{library_name_past}/details/{sample}/4.consensus/{sample}.consensus.fa"\
                .format(resultdir=library_dict["result_dir"],
                library_name_past=lib_name, 
                sample=samp)
            try:
                with open(consensus_fa, "rt") as fh:
                    next(fh)
                    gh.write(">{}_{}\n".format(samp, lib_name))
                    gh.write(fh.read())
            except Exception:
                print("WARNING - {} - no such file?".format(consensus_fa))
gh.close()
# pdb.set_trace()

# WORKENV work.env print
"""
export LIBRARY_NAME=20211231-SC2NGS
export RESULT_DIR=/home/id_seq/WORK/IDseqV2/20211231-SC2NGS
export RAWDATA_DIR=/home/id_seq/ID-seq/Share/20211231-SC2NGS-T9E9S9T9S
export HOME=/home/id_seq
"""
print("export LIBRARY_NAME={}\nexport RESULT_DIR={}\nexport RAWDATA_DIR={}\nexport HOME=/home/id_seq".\
    format(
    library_dict["library_name"],
    result_path,
    rawdata_dir
    ))
