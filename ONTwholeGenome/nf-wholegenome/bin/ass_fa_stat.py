#!/usr/bin/env python
# -*- coding: utf8 -*-

"""
Author: mengxf
Date: 2024.05.23
"""

import math
import sys
from Bio import SeqIO
import pandas as pd
from matplotlib.pyplot import MultipleLocator
from matplotlib import pyplot as plt
import click
import logging
import matplotlib

matplotlib.use('Agg')

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M', stream=sys.stdout)


class STAT:
    """
    fasta sequence 统计
    """
    def __init__(self, p_fa, prefix, outdir, asse_info):
        """
        @param p_fa:    输入fasta文件
        @param prefix:  样本前缀
        @param outdir:  输出目录
        @param asse_info:    flye 生成 assembly_info.txt 文件
        """
        self.p_fa = p_fa
        self.prefix = prefix
        self.out = outdir
        self.asse_info = asse_info
        self.Length = []

    def fa_stat(self):
        """统计质控指标"""
        global l50, l75, l90

        base_sum, gc_sum, at_sum = 0, 0, 0
        value_sum, n50, n75, n90, read_num = 0, 0, 0, 0, 0

        for record in SeqIO.parse(open(self.p_fa), "fasta"):
            num_at = record.seq.count("A") + record.seq.count("T")
            num_gc = record.seq.count("G") + record.seq.count("C")
            gc_sum += num_gc
            at_sum += num_at
            base_sum += len(record.seq)  # 计算序列总碱基数
            self.Length.append(len(record.seq))  # 序列长度列表
        self.Length.sort(reverse=True)

        n50_pos, n75_pos, n90_pos = base_sum * 0.5, base_sum * 0.75, base_sum * 0.9
        for value in self.Length:
            value_sum += value
            read_num += 1

            if n50 == 0 and n50_pos <= value_sum:
                n50 = value
                l50 = read_num
            if n75 == 0 and n75_pos <= value_sum:
                n75 = value
                l75 = read_num
            if n90 == 0 and n90_pos <= value_sum:
                n90 = value
                l90 = read_num

        # 统计长度信息
        max_length, mini_length, num = max(self.Length), min(self.Length), len(self.Length)
        av_length = base_sum / num
        gc = gc_sum / (gc_sum + at_sum) * 100
        circ_info = self.get_circ_info()

        with open(f"{self.out}/{self.prefix}_fa.stat.txt", 'w', encoding='utf-8_sig') as w:
            w.write(
                f"序列质量统计\t成环信息\tGC含量(%)\tscaffolds数量\tscaffolds总长\t最长scaffolds长度\t最短scaffolds长度"
                f"\t平均长度\tN50\tN75\tN90\tL50\tL75\tL90\n")
            w.write(f"{self.prefix}\t{circ_info}\t{gc:.2f}\t{num:,}\t{base_sum:,}\t{max_length:,}\t{mini_length:,}"
                    f"\t{av_length:.2f}\t{n50:,}\t{n75:,}\t{n90:,}\t{l50}\t{l75}\t{l90}\n")

    def get_circ_info(self):
        """统计成环信息"""
        df = pd.read_table(self.asse_info, sep='\t')
        circ_infos = []
        for idx in df.index:
            seq_name = df.loc[idx, '#seq_name']
            length = df.loc[idx, 'length']
            circ = df.loc[idx, 'circ.']
            if circ == 'Y':
                circ_string = 'circular'
            else:
                circ_string = 'linear'
            name_length_circular = f"{seq_name}:{circ_string}:{length}"
            circ_infos.append(name_length_circular)
        return ';'.join(circ_infos)

    def length_plot(self):
        """画长度分布频数直方图"""
        new_length = sorted(self.Length)
        max_length = int(max(new_length))
        logging.info(f"95%,{int(len(new_length) * 0.95)},{new_length[int(len(new_length) * 0.95)]}")

        # 数据分成40个组
        group_num = 40

        # x轴刻度间隔自适应
        x_gap = int(max_length / group_num)
        x_index = len(str(x_gap)) - 1
        x_gap = int(math.ceil(x_gap / math.pow(10, x_index)) * math.pow(10, x_index))  # 计算幂 math.pow

        seq_group = [i for i in range(0, max_length, x_gap)]
        labels = []

        # 生成x轴刻度标签：数据区间
        for i in range(0, max_length, x_gap):
            if len(labels) == 0:
                labels.append(f"<{int((i + x_gap))}")
            elif len(labels) + 2 == len(seq_group):
                labels.append(f">{int(i)}")
                break
            else:
                labels.append(f"{int(i)}-{int((i + x_gap))}")

        data = pd.DataFrame(new_length, columns=['length'])
        df = data['length'].groupby(
            pd.cut(data['length'], bins=seq_group, labels=labels)).count()  # 计数
        max_count = max(df.to_list())

        word_ticks = list(df.index)
        num_ticks = [i for i in range(len(df))]

        # 初始化一张图
        plt.figure(figsize=(18, 9))
        p1 = plt.bar(range(len(df)), df, color="#4d97cd")
        plt.xlabel('Length of scaffolds(bp)')
        plt.ylabel('Number of reads')
        plt.xticks(rotation=45, ticks=num_ticks, labels=word_ticks)  # 坐标倾斜 ,并自定义横坐标
        plt.tick_params(labelsize=6)
        plt.ylim([0, max_count * 1.1])  # 限制纵坐标
        plt.bar_label(p1, label_type='edge', size=6)  # 数据标签，edge表示将数据值标签放在柱子顶端

        # y轴刻度间隔自适应
        y_group = 5
        y_gap = int(max_count / y_group)  # 千百十位取整
        y_index = len(str(y_gap)) - 1
        y_gap = int(math.ceil(y_gap / math.pow(10, y_index)) * math.pow(10, y_index))
        if y_gap == 0:
            y_gap = 1

        logging.info(f"X轴间隔:{x_gap} ;Y轴间隔:{y_gap}")

        ax = plt.gca()
        ax.yaxis.set_major_locator(MultipleLocator(y_gap))  # y轴间隔刻度

        # 图片输出
        plt.savefig(f"{self.out}/{self.prefix}.length.png", dpi=600)  # 清晰度


# 参数
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-fa', '--fasta', required=True, type=click.Path(), help="fasta路径")
@click.option('-a', '--assembly_info', required=True, type=click.Path(), help="组装信息文件")
@click.option('-p', '--prefix', required=False, type=click.STRING, default='scaffolds', show_default=True,
              help="输出文件前缀")
@click.option('-d', '--outdir', required=False, type=click.STRING, default='.', show_default=True, help="输出文件夹")
def main(fasta, prefix, outdir, assembly_info):
    """
    flye 组装 contig 统计
    1. 统计组转序列相关质控(n50...)\n
    2. 绘制序列分布图
    """
    project = STAT(fasta, prefix, outdir, assembly_info)
    project.fa_stat()
    project.length_plot()


if __name__ == "__main__":
    main()
