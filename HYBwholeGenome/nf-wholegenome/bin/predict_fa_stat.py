#!/usr/bin/env python

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


class STAT():
    def __init__(self, p_fa, prefix, outdir):
        self.p_fa = p_fa
        self.prefix = prefix
        self.out = outdir
        self.list_contigs = []
        self.Length = []

    def fa_stat(self):

        for record in SeqIO.parse(open(self.p_fa), "fasta"):
            self.Length.append(len(record.seq))  # 序列长度列表

        self.Length.sort(reverse=True)

    def length_plot(self):
        # 画长度分布频数直方图
        new_length = sorted(self.Length)
        max_length = int(max(new_length))
        logging.info(f"95%,{int(len(new_length)*0.95)},{new_length[int(len(new_length)*0.95)]}")

        # 数据分成40个组
        group_num = 40

        # x轴刻度间隔自适应
        x_gap = int(max_length/group_num)
        x_index = len(str(x_gap))-1
        x_gap = int(math.ceil(x_gap/math.pow(10, x_index))*math.pow(10, x_index))  # 计算幂 math.pow

        seq_group = [i for i in range(0, max_length, x_gap)]
        labels = []

        # 生成x轴刻度标签：数据区间
        for i in range(0, max_length, x_gap):
            if len(labels) == 0:
                labels.append(f"<{int((i+x_gap))}")
            elif len(labels)+2 == len(seq_group):
                labels.append(f">{int(i)}")
                break
            else:
                labels.append(f"{int(i)}-{int((i+x_gap))}")

        data = pd.DataFrame(new_length, columns=['length'])
        df = data['length'].groupby(
            pd.cut(data['length'], bins=seq_group, labels=labels)).count()  # 计数
        # while df.iloc[-1] == 0:  #删除最后一行是'0'的行
        #     df.drop(df.tail(1).index,inplace=True)
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
        plt.ylim([0, max_count*1.1])  # 限制纵坐标
        plt.bar_label(p1, label_type='edge', size=6)  # 数据标签，edge表示将数据值标签放在柱子顶端

        # y轴刻度间隔自适应
        y_group = 5
        y_gap = int(max_count/y_group)   # 千百十位取整
        y_index = len(str(y_gap))-1
        y_gap = int(math.ceil(y_gap/math.pow(10, y_index))*math.pow(10, y_index))
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
@click.option('-p', '--prefix', required=False, type=click.STRING, default='scaffolds', show_default=True, help="输出文件前缀")
@click.option('-d', '--outdir', required=False, type=click.STRING, default='.', show_default=True, help="输出文件夹")
def main(fasta, prefix, outdir):
    """
    flye 组装 contig 统计
    1. 统计组转序列相关质控(N50...)\n
    2. 绘制序列分布图
    """
    project = STAT(fasta, prefix, outdir)
    project.fa_stat()
    project.length_plot()


if __name__ == "__main__":
    main()
