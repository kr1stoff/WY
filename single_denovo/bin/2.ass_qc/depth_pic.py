import sys 

if len(sys.argv) != 4:
    print('Usage: python %s <in.depth> <outdir> <reference>'%sys.argv[0])
    exit()

import os
import re
from Bio import SeqIO
from matplotlib import pyplot as plt 
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')

depthfile = sys.argv[1]
outdir = os.path.abspath(sys.argv[2])
ref = os.path.abspath(sys.argv[3])


#获取最长的基因组长度
length_list = []
for record in SeqIO.parse(open(ref), "fasta"):
    length_list.append(len(record.seq))


length = max(length_list)
length2 = length + 1 
win_num = 50

step = length2/win_num

read_len = 50 #默认read长度为75，用于估计reads数，虽然存在软切的情况会导致实际read长度小于75，计算出来的reads数会小于实际值，但不会影响太大
species_list = [[] for i in range(win_num)]
all_depth_list = []

with open(depthfile) as f:
    for line in f:
        line = line.strip()
        lines = line.split('\t')
        species, position, depth = lines[0], int(lines[1]), int(lines[2])
        all_depth_list.append(depth)

        position = position % length2
        index = int(position/step)

        species_list[index].append(depth)

out_plot_txt = '{0}/bed.depth.stat.txt'.format(outdir)
fo = open(out_plot_txt, "w")

win_dep_list = [] #用于记录每个窗口的平均深度
win_read_num_list = [] #用于记录每个窗口的总reads数
win_depth_m_list = [] #用于记录每个窗口深度中位数
for ind in range(len(species_list)):
    depth_list = species_list[ind]
    mean = np.mean(depth_list) if len(depth_list) > 0 else 0
    median = np.median(depth_list) if len(depth_list) > 0 else 0
    read_num = math.ceil(sum(depth_list)/read_len) if len(depth_list) > 0 else 0
    win_dep_list += [mean]
    win_read_num_list += [read_num]
    win_depth_m_list += [median]
    #print ("%s\t%s\t%s\t%s" %(int(ind*step),mean,median,read_num))
    fo.write("%s\t%s\t%s\t%s\n" %(int(ind*step),mean,median,read_num))

fo.close

out_stat = '{0}/bam.stat.txt'.format(outdir)
fo2 = open(out_stat, "w")

total_map_posi = len([i for i  in all_depth_list if i >=1])
total_map_posi_5 = len([i for i  in all_depth_list if i >=5])
total_map_posi_30 = len([i for i  in all_depth_list if i >=30])
total_map_posi_100 = len([i for i  in all_depth_list if i >=100])

coverage = "%.2f" % (100*total_map_posi/length)
total_map_posi_5_per = "%.2f" % (100*total_map_posi_5/total_map_posi)
total_map_posi_30_per = "%.2f" % (100*total_map_posi_30/total_map_posi)
total_map_posi_100_per = "%.2f" % (100*total_map_posi_100/total_map_posi)

total_map_base = sum(all_depth_list)  #比对上的总的碱基

avergage_depth = total_map_base/length  #平均深度
avergage_depth_20per = avergage_depth * 0.2
avergage_depth = "%.2f" % avergage_depth

total_map_posi_20_per = len([i for i  in all_depth_list if i >avergage_depth_20per])
Uniformity = "%.2f" % (100*total_map_posi_20_per/length)

fo2.write("Depth\tCoverage(%)\tDepth>=5x(%)\tDepth>=30x(%)\tDepth>=100x(%)\tUniformity\n")
fo2.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(avergage_depth,coverage,total_map_posi_5_per,total_map_posi_30_per,total_map_posi_100_per,Uniformity))
fo2.close

outpng = '{0}/bed.depth.stat.png'.format(outdir)

max_dep = max(win_dep_list)
if max_dep < 2:
    max_dep = 2
x = [i for i in range(len(win_read_num_list))]
x_step = int(len(x)/10)

plot_x = [i for i in range(0, len(win_read_num_list)+x_step, x_step)][:-1] #需要画出来的横坐标
show_x = [i for i in range(0, length, int(length/len(plot_x)+1))] #横坐标标签
plot_x += [x[-1]]
plot_x[0] = -1
show_x += [length]
if length < 10000:
    show_x = [round(i/100)*100 for i in show_x]
elif length < 1000000:
    show_x = ['{0}K'.format(round(i/1000)) for i in show_x]
elif length < 10000000:
    show_x = ['{0}M'.format(round(i/1000000, 1)) for i in show_x]
elif length < 1000000000:
    show_x = ['{0}M'.format(round(i/1000000)) for i in show_x]
else:
    show_x = ['{0}G'.format(round(i/1000000000, 1)) for i in show_x]


# 绘图
fig = plt.figure(figsize=(8,4))
ax = plt.axes()
plt.subplots_adjust(left=0.1, bottom=0.12, right=0.9, top=0.88)
ax.bar(x, win_read_num_list, color='skyblue')
#ax.bar(x, win_read_num_list, color='c')
#ax.set_title(' '.join(species.split('_')) + '\nCoverage: ' + coverage, fontsize=14) #标题
ax2 = ax.twinx() #建立次坐标轴
ax2.set_ylim(0, max_dep)
ax2.plot(x, win_dep_list, color='r', lw=0.8, label='Average Depth')
ax2.plot(x, win_depth_m_list, color='yellow', lw=1, label='Median Depth', ls='--')
#ax2.lines[1].set_linestyle('--')
#将次坐标轴设置为红色
ax2.tick_params(axis='y', colors='red')
ax2.tick_params(axis='x', colors='w')
#ax.tick_params(axis='x', colors='w')
ax2.spines['right'].set_color('red')
ax2.yaxis.label.set_color('red')
ax.spines['bottom'].set_visible(False)
ax.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.xticks(plot_x, show_x, fontsize=7)
plt.xlim(-1, win_num)
ax.set_xlabel('Genomic position', fontsize=14)
ax.set_ylabel('Mapped reads', fontsize=12)
ax2.set_ylabel('Depth', fontsize=12)
legend = plt.legend(loc=1) #图例的位置，1为右上角

plt.savefig(outpng)
plt.cla()
plt.close(fig)
