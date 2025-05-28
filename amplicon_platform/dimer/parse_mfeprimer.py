#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-01-10
import re
import pandas as pd
import sys
import os
import click
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M', stream=sys.stdout)


def find_dirst_no_space(string):
    """输出第1个不为0的位置"""
    for index, char in enumerate(string):
        if char != ' ':
            return index+1


def stat_SNP(string):
    """统计SNP数量"""
    count = 0
    for char in string:
        if char == '.':
            count += 1
    return count


def parse_group_list(group_list: list):
    """解析组列表"""
    # 解析头 Dimer 14264: S0085NT3F1_0 x S0088NT1R1_0
    res = re.split('\s+|:', group_list[0])
    idx, p1, p2 = int(res[1]), res[3], res[5]

    # 解析引物基本信息      Score: 8, Tm = 5.41 °C, Delta G = -2.06 kcal/mol
    score = int(re.findall(r"Score: (\d+)", group_list[2])[0])
    if re.findall(r"Tm = ([\d.]+) °C", group_list[2]):
        tm = float(re.findall(r"Tm = ([\d.]+) °C", group_list[2])[0])
    else:
        tm = '-'
    delta_g = float(re.findall(
        r"Delta G = ([\+\d\.-]+) kcal/mol", group_list[2])[0])

    # 引物序列/比对状态
    # print(f' {group_list[4]}',group_list[5],group_list[6])
    seq1 = str(group_list[4]).replace('\n', '')
    seq2 = str(group_list[6]).replace('\n', '')
    align_point = str(group_list[5]).replace('\n', '')

    # 计算二聚体全长 F1全长+R全长+FR比对上， 因为没有出现”: :“
    len_align = len(re.sub('\s', '', align_point))
    length = len(re.sub('\s', '', seq1)) + \
        len(re.sub('\s', '', seq2)) - len_align
    mismatch = stat_SNP(align_point)

    # seq1 seq2 比对起始/终止, 比对上起始/终止
    seq1_start, seq1_end = find_dirst_no_space(seq1), len(seq1)
    seq2_start, seq2_end = find_dirst_no_space(seq2), len(seq2)
    align_start, align_end = find_dirst_no_space(align_point), len(align_point)

    # filter
    if seq1_start == 1 and align_start == seq2_start and seq1_end == align_end:
        flag = 'Y'
    elif seq2_start == 1 and seq1_start == align_start and align_end == seq2_end:
        flag = 'Y'
    else:
        flag = 'N'

    seq1 = re.sub('\s', '', seq1)
    seq2 = re.sub('\s', '', seq2)

    # 输出
    # print(idx,f'{p1}::{p2}',score,tm,delta_g,seq1,seq2,len_align,snp,seq1_start,seq1_end,seq2_start,seq2_end,align_start,align_end,sep='\t')
    print(idx, f'{p1}::{p2}', score, tm, delta_g, seq1, seq2, len(
        seq1), len(seq2), len_align, mismatch, flag, sep='\t')


def parse_file(file):
    """遍历mfeprimer文件"""
    print("#index\tdimer\tscore\ttm\tdelta_g\tseq1\tseq2\tlen1\tlen2\tlength\tmismatch\tflag")
    with open(file, 'r', encoding='utf-8') as f:
        group_list = []
        group_flag = False
        for line in f:
            if not line:
                continue
            # 解析头 Dimer 14264: S0085NT3F1_0 x S0088NT1R1_0
            group_pattern = 'Dimer \d+: .*? x .*?'
            if re.search(group_pattern, line):
                group_flag = True
            if group_flag:
                group_list.append(line)
            if len(group_list) == 7:
                parse_group_list(group_list)

                # 重置遍历标志
                group_list = []
                group_flag = False


# =========================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--file', required=True, type=click.Path(), help="mfeprimer软件，输出txt")
def main(file):
    """
    解析 mfeprimer 软件输出的二聚体json
    [输出格式]: 序号 二聚体 得分 tm值 自由能 序列1 序列2 序列1长度 序列2长度 overlap长度 错配 是否为二聚体
    """
    parse_file(file)


if __name__ == '__main__':
    main()
