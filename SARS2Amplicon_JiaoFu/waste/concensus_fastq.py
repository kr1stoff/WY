#!/usr/bin/env python

import sys
import os
from collections import Counter


# check input
assert len(sys.argv) == 3, "ERROR - concensus_fastq.py - Need 2 Arguments, <1 Concensus FASTA> <2 Samtools mpileup>"
assert os.path.exists(sys.argv[1]), "ERROR - concensus_fastq.py - <1 Concensus FASTA> not exists!"
assert os.path.exists(sys.argv[2]), "ERROR - concensus_fastq.py - <2 Samtools mpileup> not exists!"
# assignment
concensus_fa_file   = sys.argv[1]
mpileup_file        = sys.argv[2]
output_file         = concensus_fa_file.replace("fa", "fq")


# main
mpileup_dict = dict()
with open(mpileup_file) as f:
    for line in f:
        _, pos, _, depth, _, base_qual = line.strip().split("\t")
        # gap
        if depth == "0":
            continue
        max_base_qual = str(max(map(lambda x:ord(x) - 33, base_qual)))
        mpileup_dict[pos] = max_base_qual
# most popular
counter_dict = Counter(mpileup_dict.values())
most_pop_qual = max(counter_dict, key=counter_dict.get)
# concensus
with open(concensus_fa_file) as f:
    header  = next(f)
    header  = header.strip().replace(">", "")
    context = ""
    for line in f:
        context += line.strip()
# whole bases number
whole_length = len(context)
# synthesize Q-value
qvalue_line = ""
for i in range(1, whole_length+1):
    if str(i) in mpileup_dict:
        qvalue_line += chr(int(mpileup_dict[str(i)]) + 64)
    else:
        qvalue_line += chr(int(most_pop_qual) + 64)

with open(output_file, "wt", encoding="utf-8", newline="") as g:
    # "@SARS-CoV-2/Beijing/XG0010/2022"
    g.write("@SARS-CoV-2/Beijing/{}/2022\n".format(header))
    g.write(context + "\n")
    g.write("+\n")
    g.write(qvalue_line + "\n")
