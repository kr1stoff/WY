#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import re
from argparse import ArgumentParser
from goatools import obo_parser, base
import pandas as pd


help_info = 'This script is convert the GO annotation with GO annotation with wego format, ' \
            'TSV format (interproscan.sh and eggNOG)  to GENE2TERM format'
def make_options():
    parser = ArgumentParser(description=help_info)
    parser.add_argument("-i", '--in', action='store', dest='inputfile',
                        help='path of input file',
                        default=False)
    parser.add_argument("-db", '--database', action='store', dest='db',
                        help='path of go-basic.obo',
                        default=False)

    avg = parser.parse_args()
    in_file = avg.inputfile
    db = avg.db
    return in_file,db


#
def gene_go_query(goterm_list):
    '''
    :param goterm_list: eg. goterm_list = ['GO:0042026', 'GO:0006457', 'GO:0016887', 'GO:0005524']
    :return:
    '''
    all_goterm_list = goterm_list
    for goterm in goterm_list:
        goterm_obj = obofile.query_term(goterm)
        if goterm_obj != None:
            parents_list = list(goterm_obj.get_all_parents())
            all_goterm_list = list(set(all_goterm_list + parents_list))
        else:
            all_goterm_list.remove(goterm)
            continue
    # remove top level
    # GO:0003674 molecular_function,GO:0005575 cellular_component,GO:0008150 biological_process
    all_goterm_obj_list = [obofile.query_term(go_term) for go_term in all_goterm_list if go_term != "GO:0003674" and go_term != "GO:0005575" and go_term != "GO:0008150" ]
    all_goterm_detail_list = ['\t'.join([go_term.id, "level_" + str(go_term.level), go_term.namespace, go_term.name]) for go_term in all_goterm_obj_list if go_term is not None]
    return all_goterm_detail_list
    
#解析输入文件
def parser_GO_infile(in_file):
    '''
    解析传入的emapper.annotations ，interpro结果
    eggNOG 的GO注释结果存在较多的废弃的term
    '''
    GO_dict = {}
    with open(in_file,'r',encoding='utf-8') as f:
        for line in f:
            if line.startswith("#"): continue
            if not re.findall(r'GO:\d{7}', line): continue #跳过空集
            gene_id = line.strip().split('\t')[0]
            if gene_id not in GO_dict.keys():
                GO_dict[gene_id] = []

            gene_list = re.findall(r'GO:\d{7}', line)
            tmp_list = list(set(GO_dict[gene_id] + gene_list))  #去重
            GO_dict[gene_id] = tmp_list

    # # gene_GO_dict :
    # # {'GENE1': ['GO:0003824'], 'GENE2': ['GO:0008080'],'GENE3': ['GO:0042026', 'GO:0006457', 'GO:0005524'], ... }

    out_file.write('基因ID\tGO_term\t级别\t功能分类\t功能描述\n')
    for gene_id in GO_dict.keys():
        goterm_detail_list = gene_go_query(GO_dict[gene_id])
        for GO_term in goterm_detail_list:
            go_fetch = GO_term.strip().split('\t')
            out_file.write(f"{gene_id}\t{go_fetch[0]}\t{go_fetch[1]}\t{go_fetch[2]}\t{go_fetch[3]}\n")



# 根据分类统计基因个数
def go_anno_stat():
    '''
    Summarize and plot using GO annotations at the level 2
    '''
    df_go_anno = pd.read_csv(out_file_path, sep='\t',header=0,encoding='utf-8')
    df_go_anno_level2 = df_go_anno[df_go_anno['级别'] == "level_1"]
    df_level2_counts = df_go_anno_level2.groupby(["GO_term", "功能分类", "功能描述"], as_index=False)['基因ID'].count()
    df_level2_counts = df_level2_counts.rename(columns={'基因ID': '基因数量'})
    df_level2_counts.to_csv(stat_file_path, sep='\t', index=False)



#####################################################################
# main over
if len(sys.argv) == 1:
    print(help_info)
    sys.exit()
else:
    in_file,database = make_options()
    out_dir = os.path.dirname(in_file)
    os.makedirs(os.path.join(out_dir,'GO'),exist_ok=True)

    # 加载/下载数据库
    if not database:
        database = os.path.join(out_dir,'GO','go-basic.obo')
        base.download_go_basic_obo(database)
    obofile = obo_parser.GODag(database, load_obsolete=False)
    
    #注释输出文件名    
    out_file_path = os.path.join(out_dir, 'GO','GO_anno.xls')
    stat_file_path = os.path.join(out_dir, 'GO','GO_anno_stats.xls')

    # 注释合并
    out_file = open(out_file_path, "w+")
    parser_GO_infile(in_file)
    go_anno_stat()
    out_file.close()