#!/usr/bin/env python

from operator import le
import re
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
        """统计质控指标"""
        BaseSum, GC_sum, AT_sum = 0, 0, 0
        ValueSum, N50, N75, N90, read_num = 0, 0, 0, 0, 0
        circ_infos = []

        for record in SeqIO.parse(open(self.p_fa), "fasta"):
            circ_infos.append(self.get_contig_circ_info(record.description))
            num_at = record.seq.count("A")+record.seq.count("T")
            num_gc = record.seq.count("G")+record.seq.count("C")
            GC_sum += num_gc
            AT_sum += num_at
            BaseSum += len(record.seq)  # 计算序列总碱基数
            self.Length.append(len(record.seq))  # 序列长度列表
        self.Length.sort(reverse=True)

        N50_pos, N75_pos, N90_pos = BaseSum*0.5, BaseSum*0.75, BaseSum*0.9
        for value in self.Length:
            ValueSum += value
            read_num += 1

            if N50 == 0 and N50_pos <= ValueSum:
                N50 = value
                L50 = read_num
            if N75 == 0 and N75_pos <= ValueSum:
                N75 = value
                L75 = read_num
            if N90 == 0 and N90_pos <= ValueSum:
                N90 = value
                L90 = read_num

        # 统计长度信息
        max_length, mini_legth, num = max(self.Length), min(self.Length), len(self.Length)
        av_length = BaseSum/num
        gc = GC_sum/(GC_sum+AT_sum)*100
        circ_info_string = ";".join(circ_infos)

        with open(f"{self.out}/{self.prefix}_fa.stat.txt", 'w', encoding='utf-8_sig')as w:
            w.write(
                "序列质量统计\t成环信息\tGC含量(%)\tscaffolds数量\tscaffolds总长\t最长scaffolds长度"
                "\t最短scaffolds长度\t平均长度\tN50\tN75\tN90\tL50\tL75\tL90\n")
            w.write(f"{self.prefix}\t{circ_info_string}\t{gc:.2f}\t{num:,}\t{BaseSum:,}\t{max_length:,}"
                    f"\t{mini_legth:,}\t{av_length:.2f}\t{N50:,}\t{N75:,}\t{N90:,}\t{L50}\t{L75}\t{L90}\n")

    def get_contig_circ_info(self, desc):
        """
        从给定的描述信息中提取contig的名称、长度和成环状态。

        参数:
        desc (str): 一个包含contig描述信息的字符串，预期格式为：“<名称> length=<长度> circular=<true/false>”。

        返回值:
        str: 如果符合预期格式且提取成功，则返回一个包含contig名称、成环状态（圆形或线性）和长度的字符串，
            格式为“<名称>:<成环状态>:<长度>”。如果提取失败，则返回一个空字符串，并记录一条警告日志。
        """
        pattern = re.compile("(.*?) length=(\d+).*?circular=(true|false)")
        res = re.findall(pattern, desc)  # record.description
        if res:
            name, length, circ_bool = res[0]
            circ = "circular" if circ_bool == "true" else "linear"
            return f"{name}:{circ}:{length}"
        else:
            logging.warning(
                f"record.description 不符合 (.*?) length=(\d+).*?circular=(true|false) 格式: {desc}")
            return ""

    def length_plot(self):
        """画长度分布频数直方图"""
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
    1. 统计组转序列相关质控(N50...)\n
    2. 绘制序列分布图
    """
    project = STAT(fasta, prefix, outdir)
    project.fa_stat()
    project.length_plot()


if __name__ == "__main__":
    main()
