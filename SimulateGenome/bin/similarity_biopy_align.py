#!/usr/bin/env python
#算一下病毒还行,几十K. 细菌和更大基因组算不了

import os
import sys
from Bio import SeqIO
from Bio import Align


def get_args():
    if len(sys.argv) != 3:
        print(f"Usage:\n    python {os.path.basename(sys.argv[0])} <fasta1> <fasta2>")
        sys.exit(1)
    return sys.argv[1:]

def align_score(fa1, fa2):
    aligner = Align.PairwiseAligner()
    """差异分数除fa2长度,相似度"""
    tmp1 = list(SeqIO.parse(fa1, "fasta").records)[0]
    tmp2 = list(SeqIO.parse(fa2, "fasta").records)[0]
    return aligner.score(tmp1.seq, tmp2.seq) / len(tmp2.seq)

if __name__ == "__main__":
    fa1, fa2 = get_args()
    print(f"输入的序列相似度为：{align_score(fa1, fa2)}") #:.2%
