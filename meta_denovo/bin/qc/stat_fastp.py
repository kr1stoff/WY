#!/usr/bin/python
#-*- coding:utf-8 -*-
#Author: wjchen@visionmedicals.com
#Creat time: 2019-09-20 18:08:47
#Last modified: 2020-07-06 21:49:28

import os
import sys
import re
import json

if len(sys.argv) != 3:
    print('Usage: python %s <in.json> <out.xls>'%sys.argv[0])
    sys.exit(-1)


infile, outfile = sys.argv[1:]

sample = os.path.basename(infile).split('.')[0]

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
OUT = open(outfile, 'w')
OUT.write('sample\traw_reads\tclean_reads\tclean_rate\tadapter_rate\tduplication_rate\traw_Q20\traw_Q30\traw_GC\tclean_Q20\tclean_Q30\tclean_GC\n')
OUT.write('{0}\t{1}\t{2}\t{3:.1f}\t{4:.1f}\t{5:.1f}\t{6:.1f}\t{7:.1f}\t{8:.1f}\t{9:.1f}\t{10:.1f}\t{11:.1f}\n'.format(sample, raw_reads, clead_reads, clead_rate, adapter_rate, duplication, raw_q20, raw_q30, raw_gc_content, clean_q20, clean_q30, clean_gc_content))

