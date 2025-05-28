#!/usr/bin/python
#-*- coding:utf-8 -*-
#Author: wjchen@visionmedicals.com
#Creat time: 2019-09-20 18:08:47
#Last modified: 2020-07-06 21:49:28

import os
import sys
import re
import json

if len(sys.argv) != 4:
    print('Usage: python %s <in.json> <out.xls> <out.png>'%sys.argv[0])
    sys.exit(-1)

#import matplotlib.pyplot as plt
#import matplotlib as mpl 
 
  
#import seaborn as sns  
#sns.set_style('darkgrid')

infile, outfile, outpng = sys.argv[1:]

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
OUT.write('#sample\traw_reads\tclean_reads\tclean_rate\tadapter_rate\tduplication_rate\traw_Q20\traw_Q30\traw_GC\tclean_Q20\tclean_Q30\tclean_GC\n')
OUT.write('{0}\t{1}\t{2}\t{3:.1f}\t{4:.1f}\t{5:.1f}\t{6:.1f}\t{7:.1f}\t{8:.1f}\t{9:.1f}\t{10:.1f}\t{11:.1f}\n'.format(sample, raw_reads, clead_reads, clead_rate, adapter_rate, duplication, raw_q20, raw_q30, raw_gc_content, clean_q20, clean_q30, clean_gc_content))

#绘制碱基质量图
#A = data["read1_before_filtering"]["quality_curves"]["A"]
#T = data["read1_before_filtering"]["quality_curves"]["T"]
#C = data["read1_before_filtering"]["quality_curves"]["C"]
#G = data["read1_before_filtering"]["quality_curves"]["G"]

#x = [i for i in range(len(A))]


#A = [float(i)+2 for i in A]
#T = [float(i)+2 for i in T]
#C = [float(i)+2 for i in C]
#G = [float(i)+2 for i in G]
#win_num = 75
#win_step = len(A)/win_num
#plot_x = [i + win_step/2 for i in [b*win_step for b in range(9)]]
#plot_x = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]
#plot_x += [i + win_step/2 for i in [b*win_step for b in range(0, win_num+1, 2)][5:]]
#show_x = [int(i/win_step)+1 for i in plot_x]
#x = [i+win_step/2 for i in range(len(A))]
#fig = plt.figure(figsize=(11.25,6))
#ax = plt.axes()
#plt.plot(x, A, label='A', lw=1)
#plt.plot(x, T, label='T', lw=1)
#plt.plot(x, C, label='C', lw=1)
#plt.plot(x, G, label='G', lw=1)
#plt.ylim(0, 36)
#plt.xlim(0, len(A))
#plt.title('Quality scores cross all bases (Sanger / Illumina 1.9 encoding)')
#plt.xlabel('Position in read (bp)')
#plt.xticks(plot_x, show_x)
#plt.yticks([i for i in range(0, 37, 2)])
#plt.axhspan(0, 20, facecolor='red', alpha=0.3)
#plt.axhspan(20, 30, facecolor='wheat', alpha=0.8)
#plt.axhspan(30, 36, facecolor='limegreen', alpha=0.4)
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.tick_params(axis="x", bottom=False)
#ax.tick_params(axis="y", left=False)
#for k in range(0, win_num, 2):
#    plt.axvspan(win_step*k, win_step*(k+1), facecolor='white', alpha=0.3)
#plt.legend(loc=4)
#plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=0.1, wspace=0.1)
#plt.subplots_adjust(left=0.03, bottom=0.08, right=0.97, top=0.92)
#plt.savefig(outpng)
