#!/home/earthtest/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/8/11 18:00

import os
import sys
import pandas as pd 
from Bio import SeqIO


if len(sys.argv) != 5:
    print(f"Useage:\n\tpython {sys.argv[0]} <depth_file> <fa> <out> <per num>\n")
    print("Description:\n\t根据深度文件，设置每间隔 (num) bp ,计算位点范围内的 av_depth 与 GC% ，仅限比对序列只有1条")
    exit(-1)

def main(dfile,fa,out,num):
    df=pd.DataFrame(pd.read_csv(dfile,sep='\t',encoding='utf-8_sig',names=['ref','pos','depth']))
    OUT=open(out,'w',encoding='utf-8_sig')
    OUT.write(f"ID\tStart\tEnd\tDepth\tGC\n")

    for tmp_record in SeqIO.parse(fa, 'fasta'):
        if df.iloc[0, 0]== tmp_record.id:
            record=tmp_record

    tmp_i=0
    for i in range(0,len(df),100):
        if i ==0 :continue

        av_dep=df.iloc[tmp_i:i,2].mean()
        num_gc=f"{(record.seq[tmp_i:i].count('G')+record.seq[tmp_i:i].count('C'))/int(num)*100:.2f}"

        OUT.write(f"{record.id}\t{tmp_i+1}\t{i}\t{av_dep}\t{num_gc}\n")
        tmp_i=i

    OUT.close()



if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
