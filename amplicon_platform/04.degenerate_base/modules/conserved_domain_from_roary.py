import re
from BCBio import GFF


def get_locus_locations_dict(ref_gff_file, thres_length):
    """
    参考基因组gff文件, 提取locus_tag位置信息. roary使用locus_tag代替基因
    params:
    ref_gff_file    ::  参考基因组gff文件
    thres_length    ::  locus_tag长度阈值, 默认100
    return:
    dict_locus_location ::  {'GCF_001160345_00001': ['NZ_CQAE01000001.1', '38', '779'], ...}
    """
    dict_locus_location = {}
    for rec in GFF.parse(open(ref_gff_file)):
        if rec.features == []: continue #没有feature的chrom
        for feat in rec.features:
            quali = feat.qualifiers
            if ('locus_tag' not in quali) or ((feat.location.end-feat.location.start)<=thres_length): continue
            locus_tag = quali['locus_tag'][0]
            chrom = rec.id
            start, end = str(feat.location.start), str(feat.location.end)
            dict_locus_location[locus_tag] = [chrom, start, end]
    return dict_locus_location


def get_core_cluster_presence_dict(gene_presence_absence, ref_assembly, thres_presence):
    """
    筛选包容性95以上的核心蛋白群(clustered_proteins), 蛋白群包含各基因组gff的locus_tag
    params:
    gene_presence_absence   ::  roary结果,蛋白群出现/缺席二进制矩阵表
    ref_assembly            ::  参考/代表性序列assembly_accession
    thres_presence          ::  出勤率阈值, 默认0.95
    return:
    dict_core_cluster_presence  ::  {'kdsB': '99.78%', ...}
    """
    dict_core_cluster_presence = {}
    with open(gene_presence_absence) as f:
        #表头.有效基因组数量,参考/代表基因组所在列
        header = next(f)
        headers = header.strip().split('\t')
        valid_genome_number = len(headers) - 1
        ref_column = headers.index(ref_assembly)
        for line in f:
            llst = line.strip().split('\t')
            #[221209]gene_presence_absence.Rtab基因控空格有空格用双引号"",clustered_proteins不用, 这里不一致
            cluster = llst[0].strip('"')
            binaries = list(map(int, llst[1:]))
            bsum = sum(binaries)
            presence_percent = bsum / valid_genome_number
            if (binaries[ref_column-1] == 1) and (presence_percent >= thres_presence): #参考/代表基因组有该基因,包容性95以上
                dict_core_cluster_presence[cluster] = f'{presence_percent:.2%}'
    return dict_core_cluster_presence


def get_cluster_locus_dict(clustered_proteins, dict_core_cluster_presence, ref_assembly):
    """生成cluster和locus对照字典
    clustered_proteins              ::  roary结果,cluster和所有基因组的locus的表格
    dict_core_cluster_presence      ::  get_core_cluster输出的核心蛋白群
    ref_assembly                    ::  参考/代表性序列assembly_accession
    return:
    dict_cluster_locus  ::  {'pytH': 'GCF_001160345_03107', ...}
    """
    dict_cluster_locus = {}
    with open(clustered_proteins) as f:
        for line in f:
            cluster, loca = line.strip().split(':')
            if cluster in dict_core_cluster_presence:
                locus = re.findall(f'({ref_assembly}_\d+)', loca)[0]
                dict_cluster_locus[cluster] = locus
    return dict_cluster_locus


def create_cluster_location_presence_table(
    file_core_cluster_location_presence, 
    bed_core_cluster,
    dict_cluster_locus:dict,
    dict_locus_location:dict, 
    dict_core_cluster_presence:dict):
    """合并三个字典, 输出结果表格"""
    with open(file_core_cluster_location_presence, 'wt', encoding='utf-8', newline='') as g, \
    open (bed_core_cluster, 'wt', encoding='utf-8', newline='') as h:
        g.write('#cluster\tlocus\tchrom\tstart\tend\tpresence\n')
        for cluster, presence in dict_core_cluster_presence.items():
            locus = dict_cluster_locus[cluster]
            if locus not in dict_locus_location: continue #[221206]长度阈值过滤掉的locus_tag
            chrom, start, end = dict_locus_location[locus]
            g.write('\t'.join([cluster, locus, chrom, start, end, presence]) + '\n')
            h.write(f'{chrom}\t{start}\t{end}\t{locus}\n') #221206 BCBio.GFF use python 0-base
