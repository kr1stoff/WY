#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##导入模块，初始传递命令、变量等
import argparse
import re

parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 base_stat.py -g [genome.fasta] -d [samtools.depth] -s [depth_base.stat] -l [split_length]\n\n输入基因组 fasta 文件，以及由 Samtools 得到的基因序列各位点测序深度统计文件\n输出按基因组一定大小滑窗进行统计后的包含滑窗 reads 覆盖深度及碱基含量百分比的表格')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-g', '--genome', metavar = '[genome.fasta]', help = '输入文件，基因组 fasta 文件', required = True)
required.add_argument('-d', '--depth', metavar = '[samtools.depth]', help = '输入文件，原始位点测序深度统计文件（由 Samtools 得到）', required = True)
required.add_argument('-s', '--stat', metavar = '[depth_base.stat]', help = '输出文件，按基因组一定大小滑窗进行统计后的，包含滑窗 reads 覆盖深度及碱基含量百分比的表格', required = True)
optional.add_argument('-l', '--split', metavar = '[split_length]', type = int, default = 2000, help = '序列统计间隔，默认 2000bp 为一统计滑窗', required = False)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

##读入并处理文件
#读取基因组，并统计基因组 A、T、G、C、GC 比例

#传入的fa，记录每条序列的序列信息到 字典——genome_dict
#传入的fa，记录每条序列的长度到 字典——genome_size
genome_dict = {}
genome_size = {}
with open(args.genome, 'r') as genome_fasta:
	for line in genome_fasta:
		line = str(line.strip())
		
		if line.startswith(">"):#头
			seq = line.split('>')[1]
			genome_dict[seq] = ''
			genome_size[seq] = 0
		else:#序列
			genome_dict[seq] = genome_dict[seq] + line
			genome_size[seq] = genome_size[seq] + len(line)


result = {}
n = 0
for seq_ID,seq_seq in genome_dict.items():	#遍历序列信息
	seq_length = len(seq_seq)
	seq_start = list(range(0, seq_length, args.split))			#生成滑窗列表 0，6500，2000
	seq_end = list(range(args.split, seq_length, args.split))	# 2000,6500,2000
	seq_end.append(seq_length)


	for i in range(0, len(seq_start)):
		n = n + 1
		split_seq = seq_seq[seq_start[i]:seq_end[i]]	#按滑窗位置，切割fa
		split_len = seq_end[i] - seq_start[i]
		result[n] = [f'{seq_ID}\t{seq_start[i] + 1}\t{seq_end[i]}']
		result[n].append('0')
		
		#计算各碱基含量
		result[n].append(round(100 * len(re.findall('[GCgc]', split_seq)) / split_len, 2))
		result[n].append(round(100 * len(re.findall('[Aa]', split_seq)) / split_len, 2))
		result[n].append(round(100 * len(re.findall('[Tt]', split_seq)) / split_len, 2))
		result[n].append(round(100 * len(re.findall('[Gg]', split_seq)) / split_len, 2))
		result[n].append(round(100 * len(re.findall('[Cc]', split_seq)) / split_len, 2))

##读入 samtools 位点深度统计文件，并按滑窗统计 reads 覆盖深度
n = 0
with open(args.depth, 'r') as samtools_depth:
	line = samtools_depth.readline().strip().split('\t')
	seq_ID = line[0]
	split_start = 1
	split_end = args.split
	
	if int(line[1]) > args.split:
		for gap in range(split_start, int(line[1]) - args.split, args.split):
			n = n + 1
			result[n][1] = 0

			split_start = split_start + args.split
			split_end = split_end + args.split

	seq_end = int(line[1])
	depth = int(line[2])
	
	#遍历深度文件
	for line in samtools_depth:
		line = line.strip().split()
		if line[0] == seq_ID:
			if int(line[1]) > args.split:
				for gap in range(split_start, int(line[1]) - args.split, args.split):
					n = n + 1
					
					result[n][1] = int(depth / args.split)
					split_start = split_start + args.split
					split_end = split_end + args.split
					depth = 0
			depth = depth + int(line[2])
			seq_end = int(line[1])
		
		else:
			if seq_end == genome_size[seq_ID]:
				n = n + 1
				result[n][1] = int(depth / (seq_end - split_start + 1))
			elif genome_size[seq_ID] - split_start + 1 <= args.split:
				seq_end = genome_size[seq_ID]
				n = n + 1
				result[n][1] = int(depth / (seq_end - split_start + 1))
			else:
				n = n + 1
				result[n][1] = int(depth / args.split)
				split_start = split_start + args.split
				for gap in range(split_start, genome_size[seq_ID]- args.split, args.split):
					n = n + 1
					result[n][1] = 0
					split_start = split_start + args.split
					split_end = split_end + args.split
				n = n + 1
				result[n][1] = 0
			
			seq_ID = line[0]
			split_start = 1
			split_end = args.split
			if int(line[1]) > args.split:
				for gap in range(split_start, int(line[1]) - args.split, args.split):
					n = n + 1
					result[n][1] = 0
					split_start = split_start + args.split
					split_end = split_end + args.split
			depth = int(line[2])
			seq_end = int(line[1])
	
	if seq_end == genome_size[seq_ID]:
		n = n + 1
		result[n][1] = int(depth / (seq_end - split_start + 1))
	elif genome_size[seq_ID] - split_start + 1 <= args.split:
		seq_end = genome_size[seq_ID]
		n = n + 1
		result[n][1] = int(depth / (seq_end - split_start + 1))
	else:
		n = n + 1
		result[n][1] = int(depth / args.split)
		split_start = split_start + args.split
		for gap in range(split_start, genome_size[seq_ID]- args.split, args.split):
			n = n + 1
			result[n][1] = 0
			split_start = split_start + args.split
			split_end = split_end + args.split
		n = n + 1
		result[n][1] = 0

#输出结果
depth_base_stat = open(args.stat, 'w')
print('seq_ID\tseq_start\tseq_end\tdepth\tGC\tA\tT\tG\tC', file = depth_base_stat)
for key in range(1, len(result) + 1):
	result[key] = [str(i) for i in result[key]]
	print('\t'.join(result[key]), file = depth_base_stat)

depth_base_stat.close()
genome_fasta.close()
samtools_depth.close()
