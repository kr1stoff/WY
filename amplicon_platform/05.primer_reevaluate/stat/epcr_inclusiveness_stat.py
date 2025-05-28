#!/usr/bin/python3
# @auther: zhuzi
# @date: 2023-07-11
import sys
import click
import logging
from operator import itemgetter  # 指定列表多个索引直接提取
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M', stream=sys.stdout)


def file_iterator(file):
    """返回一个解析文本迭代器,直接迭代每1行制表符拆分的列表"""
    with open(file, 'r', encoding='utf-8') as f:
        for line in f:
            if not line:
                continue
            arg = line.strip().split('\t')
            yield arg


def parse_chr_file(file):
    """解析assembly_chromosome"""
    gene_set = set()
    args = file_iterator(file)
    for arg in args:
        gene_set.add(arg[1])  # 添加染色体名称
    return gene_set


def parse_probe_epcr(file):
    """
    解析探针epcr结果
    [Format]: { GCF: { chr: [[start,end] ,[...]] } }
    """
    probe_dic = {}
    args = file_iterator(file)
    indices = [1, 8, 9, 18]
    for arg in args:
        seqid, start1, end1, chr = itemgetter(*indices)(arg)[0:]
        start = start1 if int(start1) < int(end1) else end1  # 防止负链，位置反向
        end = end1 if int(start1) < int(end1) else start1
        tmp = [start, end]
        probe_dic.setdefault(chr, {})
        probe_dic[chr].setdefault(seqid, [])
        probe_dic[chr][seqid].append(tmp)
    return probe_dic


def parse_primer_epcr(file):
    """解析引物epcr结果"""
    primer_dic = {}
    args = file_iterator(file)
    indices = [1, 5, 6, 12, 13, 20]
    for arg in args:
        seqid, start1, end1, start2, end2, chr = itemgetter(*indices)(arg)[0:]
        start = start1 if int(start1) < int(end1) else end1  # 防止负链，位置反向
        end = start2 if int(start2) > int(end2) else end2
        tmp = [start, end]
        primer_dic.setdefault(chr, {})
        primer_dic[chr].setdefault(seqid, [])
        primer_dic[chr][seqid].append(tmp)
    return primer_dic


def judge_region(region1, region2, dic: dict, chr):
    """ 用于判断内引物是否有被外引物包含 """
    res_list = dic.get(chr, None)
    region_list = []
    if not res_list:
        return region_list
    else:
        for res in res_list:
            if int(region1) <= int(res[0]) and int(region2) >= int(res[1]):
                add_region = [int(res[0]), int(res[1])]
                region_list.append(add_region)  # 添加多个被包含的区域
    return region_list


# =========================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--id', required=True, type=click.STRING, default='test', help="菌种名称，用于做输出结果表头,默认[test]")
@click.option('-1', '--inner', required=True, type=click.Path(), help="内引物ePCR结果")
@click.option('-2', '--outside', required=True, type=click.Path(), help="内引物ePCR结果")
@click.option('-3', '--probe', required=True, type=click.Path(), help="内引物ePCR结果")
@click.option('-c', '--chr', required=True, type=click.Path(), help="assembly_chromosome")
def main(id, inner, outside, probe, chr):
    """
    统计引物设计结果每对引物的包容性\n
    判断条件：只要epcr结果中，某个GCA/F基因组中任何1条chr，满足有外引物区域包含其中1对内引物区域，且该对内引物也包含其中1条探针区域，则该基因组符合包容性条件
    """
    gene_set = parse_chr_file(chr)
    probe_hash = parse_probe_epcr(probe)
    inprimer_hash = parse_primer_epcr(inner)
    outprimer_hash = parse_primer_epcr(outside)

    len_genome = 0                # 比对上的基因组数量
    inclusive_genome_set = set()  # 包容性基因组
    for gene in gene_set:         # 遍历基因组集合
        gene_flag = False
        in_dic = inprimer_hash.get(gene, None)
        probe_dic = probe_hash.get(gene, None)
        out_dic = outprimer_hash.get(gene, None)

        if not out_dic:
            continue     # 外引物没有比对到该染色体的记录
        if not in_dic:
            continue      # 内引物没有比对到该染色体的记录
        if not probe_dic:
            continue   # 探针没有比对到该染色体的记录

        # 判断外引物是否包含内引物
        for chr in out_dic.keys():
            inregion_list = []
            for record in out_dic.get(chr, None):
                if not record:
                    continue     # 外引物没有比对到该序列的记录
                out_start, out_end = record[0:]
                region_list = judge_region(out_start, out_end, in_dic, chr)
                inregion_list = inregion_list + region_list

            # 判断被包含的内引物区域是否包含探针
            if not inregion_list:
                continue
            probe_region_list = []
            for record in inregion_list:
                in_start, in_end = record[0:]
                region_list = judge_region(in_start, in_end, probe_dic, chr)
                probe_region_list = probe_region_list + region_list

            # 包容性基因组计算
            if not gene_flag and inregion_list and probe_region_list:
                len_genome += 1
                gene_flag = True
                inclusive_genome_set.add(gene)

    # 包容性基因组统计结果打印
    genome_str = ','.join(inclusive_genome_set)
    diff_genome = ','.join(gene_set-inclusive_genome_set)  # 差集
    inclusive_rate = format(len(inclusive_genome_set)/len(gene_set)*100, '.2f')
    print('#spe_id\tlen_inclusive_genome\ttotal_genome\tinclusive_rate\tdiff_genome')
    print(id, len_genome, len(gene_set), inclusive_rate, diff_genome, sep='\t')


if __name__ == '__main__':
    main()
