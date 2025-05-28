#!/usr/bin/python
#-*- coding:utf-8 -*-
#Author: wjchen@visionmedicals.com
#Creat time: 2019-09-20 18:08:47
#Last modified: 2020-07-06 21:49:28

from encodings import utf_8_sig
import os
import sys 
import re
import json

if len(sys.argv) != 4:
    print('Usage: python %s <sample> <fq_json> <out.xls>'%sys.argv[0])
    sys.exit(-1)


sample, infile,outfile = sys.argv[1:]
#print(len(sys.argv))





with open(infile) as f:
    data = json.load(f)

raw_reads = int(data['summary']['before_filtering']['total_reads'])
clead_reads = int(data['summary']['after_filtering']['total_reads'])
clead_rate = 100*(float(clead_reads/raw_reads)) if raw_reads > 0 else 0

raw_q20 = 100*float(data['summary']['before_filtering']['q20_rate'])
raw_q30 = 100*float(data['summary']['before_filtering']['q30_rate'])
clean_q20 = 100*float(data['summary']['after_filtering']['q20_rate'])
clean_q30 = 100*float(data['summary']['after_filtering']['q30_rate'])

raw_gc_content = 100*float(data['summary']['before_filtering']['gc_content'])
clean_gc_content = 100*float(data['summary']['after_filtering']['gc_content'])

duplication = 100*float(data['duplication']['rate'])
try:
    adapter_reads = float(data['adapter_cutting']['adapter_trimmed_reads'])
    adapter_rate = 100*(adapter_reads/int(raw_reads)) if raw_reads > 0 else 0
except:
    adapter_reads = 0
    adapter_rate = 0

#输出统计结果
OUT = open(outfile,encoding='utf_8_sig',mode= 'w')
OUT.write('样本编号\t原始总reads数\t过滤后reads数\t过滤比例\t接头比例\t重复reads比例\t原始数据Q20\t原始数据Q30\t原始数据GC%\t过滤数据Q20\t过滤数据Q30\t过滤数据GC%\n')
OUT.write('{0}\t{1}\t{2}\t{3:.1f}\t{4:.1f}\t{5:.1f}\t{6:.1f}\t{7:.1f}\t{8:.1f}\t{9:.1f}\t{10:.1f}\t{11:.1f}\n'.format(sample, raw_reads, clead_reads, clead_rate, adapter_rate, duplication, raw_q20, raw_q30, raw_gc_content, clean_q20, clean_q30, clean_gc_content))

