#!/usr/bin/env python
################################################################################
# prgram:   vcf_filter_4nextclade.py
# data:     2022/02/17
# input:    sample.vcf
# output:   sample.nextclade.vcf
################################################################################
import sys
import os
import pdb


# function
def verify_flank(position:int, step:int=10):
    boolean = False
    for pos in nextclade_pos_list:
        if (pos - step) <= position <= (pos + step):
            boolean = True
            break
    return boolean

# check input
assert len(sys.argv) == 2, "ERROR - vcf_filter_4nextclade.py - Need 1 Arguments, <VCF File>"
assert os.path.exists(sys.argv[1]), "ERROR - vcf_filter_4nextclade.py - <1 Concensus FASTA> not exists!"
# assignment IO
invcf = sys.argv[1]
output_file = invcf.replace(".vcf", ".nextclade.vcf")

# configure file position
program_dir         = sys.path[0]
nextclade_pos_file  = os.path.split(program_dir)[0] + os.sep + "libs/nextclade_all.pos"
assert os.path.isfile(nextclade_pos_file), "ERROR - vcf_filter_4nextclade.py - Nextclade position file not exists!"

# nextclade 
nextclade_pos_list = list()
with open(nextclade_pos_file) as f:
    for line in f:
        pos = int(line.strip())
        if pos not in nextclade_pos_list:
            nextclade_pos_list.append(pos)

# main
#`CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  XG0008
with open(invcf) as fh, open(output_file, "wt", encoding="utf-8",newline="") as gh:
    for line in fh:
        if line.startswith("#"):
            gh.write(line)
        else:
            linelist = line.strip().split()
            # verify nexclade position flank
            if not verify_flank(int(linelist[1]), step=len(linelist[3])):
                continue
            else:
                gh.write(line)
