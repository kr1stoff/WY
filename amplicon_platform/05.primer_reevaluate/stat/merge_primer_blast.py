#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import click
from itertools import product


def parse_table(file):
    dic = {}
    ref_set = set()
    with open(file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if not line:
                continue
            # name,ref,taxid,fr,strand1,start1,end1,mismatch1,gap1,start_gap1,fr,strand2,start2,end2,mismatch2,gap2,start_gap2,all_len,total_len,strand,same_diff
            # qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus staxid sscinames scomnames
            args = line.strip().split('\t')
            dic.setdefault(args[1], [])
            dic[args[1]].append(args)
            ref_set.add(args[1])
    return dic, ref_set


def judge_region(list1: list, list2: list, list3=None):
    """判断区间是否包含"""
    flag = False
    if list3:
        # 筛选出内引物区间完全 被外引物区间 包含的情况
        if int(list1[5]) >= int(list2[5]) and int(list1[13]) <= int(list2[13]):
            if int(list3[8]) < int(list3[9]):   # blast结果负链的话，start 和 end 的位置会换过来
                # 筛选出探针区间完全 被内引物区间 包含的情况
                if int(list3[8]) >= int(list1[5]) and int(list3[9]) <= int(list1[13]):
                    flag = True

            elif int(list3[8]) > int(list3[9]):
                if int(list3[9]) >= int(list1[5]) and int(list3[8]) <= int(list1[13]):
                    flag = True

    else:
        # 筛选出内引物区间完全 被外引物区间 包含的情况
        if int(list1[5]) >= int(list2[5]) and int(list1[13]) <= int(list2[13]):
            flag = True

    return flag


# === Option ===========================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-1', '--inprimer_blast', required=True, type=click.Path(exists=True), help="内引物 inprimer_blast.same 结果")
@click.option('-2', '--outprimer_blast', required=True, type=click.Path(exists=True), help="外引物 outprimer_blast.same 结果")
@click.option('-3', '--probe_blast', required=False, type=click.Path(exists=True), help="探针 probe_blast.same 结果(可选项)")
def main(inprimer_blast, outprimer_blast, probe_blast):
    in_dic, in_set = parse_table(inprimer_blast)
    out_dic, out_set = parse_table(outprimer_blast)
    if probe_blast:
        probe_dic, probe_set = parse_table(probe_blast)

        re_header = ['#name', 'id', 'taxid', 'in_start', 'in_end', 'in_all_len', 'in_total_len', 'in_strand', 'in_SNP1', 'in_SNP2', 'out_start', 'out_end', 'out_all_len', 'out_total_len',
                     'out_strand', 'out_SNP1', 'out_SNP2', 'probe_start', 'probe_end', 'probe_length', 'probe_strand', 'probe_SNP', 'in_same_diff', 'out_same_diff', 'probe_same_diff']
        print('\t'.join(re_header))

        # 参考 id 交集
        ref_set = in_set & out_set & probe_set
        for ref in ref_set:
            in_res = in_dic.get(ref)
            out_res = out_dic.get(ref)
            probe_res = probe_dic.get(ref)

            # 排列组合
            loop_val = [in_res, out_res, probe_res]
            for i in product(*loop_val):
                li1, li2, li3 = i[0], i[1], i[2]
                if judge_region(li1, li2, li3):
                    in_SNP1 = int(li1[7]) + int(li1[8]) + int(li1[9])
                    in_SNP2 = int(li1[14]) + int(li1[15]) + int(li1[16])
                    out_SNP1 = int(li2[7]) + int(li2[8]) + int(li2[9])
                    out_SNP2 = int(li2[14]) + int(li2[15]) + int(li2[16])
                    probe_SNP = int(li3[4]) + int(li3[5])
                    print(li1[0], li1[1], li1[2], li1[5], li1[13], li1[17], li1[18], li1[19], in_SNP1, in_SNP2, li2[5], li2[13], li2[17],
                          li2[18], li2[19], out_SNP1, out_SNP2, li3[8], li3[9], li3[12], li3[17], probe_SNP, li1[20], li2[20], li3[18], sep='\t')

    else:   # 不输入探针表格
        re_header = ['#name', 'id', 'taxid', 'in_start', 'in_end', 'in_all_len', 'in_total_len', 'in_strand', 'in_SNP1', 'in_SNP2',
                     'out_start', 'out_end', 'out_all_len', 'out_total_len', 'out_strand', 'out_SNP1', 'out_SNP2', 'in_same_diff', 'out_same_diff']
        print('\t'.join(re_header))

        # 参考 id 交集
        ref_set = in_set & out_set
        for ref in ref_set:
            in_res = in_dic.get(ref)
            out_res = out_dic.get(ref)

            # 排列组合
            loop_val = [in_res, out_res]
            for i in product(*loop_val):
                li1, li2 = i[0], i[1]
                if judge_region(li1, li2):
                    in_SNP1 = int(li1[7]) + int(li1[8]) + int(li1[9])
                    in_SNP2 = int(li1[14]) + int(li1[15]) + int(li1[16])
                    out_SNP1 = int(li2[7]) + int(li2[8]) + int(li2[9])
                    out_SNP2 = int(li2[14]) + int(li2[15]) + int(li2[16])
                    print(li1[0], li1[1], li1[2], li1[5], li1[13], li1[17], li1[18], li1[19], in_SNP1, in_SNP2,
                          li2[5], li2[13], li2[17], li2[18], li2[19], out_SNP1, out_SNP2, li1[20], li2[20], sep='\t')


if __name__ == '__main__':
    main()
