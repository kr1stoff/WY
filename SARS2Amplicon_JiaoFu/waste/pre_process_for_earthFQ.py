#!/usr/bin/env python

##################################################################
## Example:
# python3 ${SCRIPTS_PATH}/utils/pre_process.py library.json
##################################################################

import sys
import os
import glob
import pdb


_self_path = os.path.basename(sys.argv[0])
# check
assert len(sys.argv) == 3, "ERROR - {} - Need 2 Arguments! <task_sample.csv> <result_dir>".format(_self_path)
assert os.path.isfile(sys.argv[1]), "ERROR - {} - Need <task_sample.csv>, e.g.S1\\tS1R1.fq,S1R2.fq".format(_self_path)
# assign arguments
task_sample = sys.argv[1]
result_dir = sys.argv[2]

# intermedia directory & sample_path_table.txt
intermedia_path = result_dir + os.sep + "intermedia"
if not os.path.isdir(intermedia_path):
    os.makedirs(intermedia_path)
sample_path_table_file  = intermedia_path + os.sep + "sample_path_table.txt"

# IO
with open(task_sample, "rt") as fh, open(sample_path_table_file, "wt", encoding="utf-8", newline="") as gh:
    for line in fh:
        llist = line.strip().split("\t")
        # check
        assert len(llist) == 3, "ERROR - {} - Check task_sample.csv, need PE FASTQ data!".format(_self_path)
        assert os.path.isfile(llist[1]), "ERROR - {} - Raw FASTQ: {} not exists!".format(_self_path, llist[1])
        assert os.path.isfile(llist[2]), "ERROR - {} - Raw FASTQ: {} not exists!".format(_self_path, llist[2])
        name, read1, read2 = llist
        sample_path_line = "{}\t{},{}\n".format(name, read1, read2)
        gh.write(sample_path_line)
