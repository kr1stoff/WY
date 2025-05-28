#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import click
from Bio import SeqIO


# COG的颜色配置
info_cog = {'A': "#FF9900\tRNA processing and modification",
            'B': "#FFB38C\tChromatin structure and dynamics",
            'C': "#FFFF00\tEnergy production and conversion",
            'D': "#CC0033\tCell cycle control, cell division, chromosome partitioning",
            'E': "#99CC33\tAmino acid transport and metabolism",
            'F': "#0000FF\tNucleotide transport and metabolism",
            'G': "#99CCFF\tCarbohydrate transport and metabolism",
            'H': "#CCFF99\tCoenzyme transport and metabolism",
            'I': "#00FF00\tLipid transport and metabolism",
            'J': "#FF0000\tTranslation, ribosomal structure and biogenesis",
            'K': "#FFCB2F\tTranscription",
            'L': "#FF9999\tReplication, recombination and repair",
            'M': "#33CC33\tCell wall/membrane/envelope biogenesis",
            'N': "#CCFFFF\tCell motility",
            'O': "#99CC00\tPosttranslational modification, protein turnover, chaperones",
            'P': "#663300\tInorganic ion transport and metabolism",
            'Q': "#003366\tSecondary metabolites biosynthesis, transport and catabolism",
            'R': "#999999\tGeneral function prediction only",
            'S': "#333333\tFunction unknown",
            'T': "#FF99CC\tSignal transduction mechanisms",
            'U': "#CCCC66\tIntracellular trafficking, secretion, and vesicular transport",
            'V': "#996699\tDefense mechanisms",
            'W': "#993399\tExtracellular structures",
            'Y': "#993366\tNuclear structure",
            'Z': "#00FFFF\tCytoskeleton"}


#解析eggnog结果
def parse_eggnog(fname):
    """
    #列出前5行eggnog结果
    Parse the eggNOG file
    1## Wed Feb 23 15:53:03 2022
    2## emapper-2.1.6
    3## emapper.py -m diamond --cpu 24 -d bact --override 
    4##
    5##query  seed_ortholog   evalue  score   eggNOG_OGs
    6##   正式结果  
    """

    res = {}
    with open(fname, 'r') as IN:
        #跳过前面4行注释行
        next(IN)
        next(IN)
        next(IN)
        next(IN)

        #第五行，为首字符带"#"的标题行
        header = next(IN).lstrip("#").strip().split("\t")
        for line in IN:
            #最后的尾行也有"#"的注释行
            if line.startswith("##"):
                pass
            
            else:
                arr = line.strip().split('\t')
                res[arr[0]] = arr[1:]
    #返回一个标题，和一个字典 key为序列名，
    return header, res


#封装程序  main
@click.command(context_settings= dict(help_option_names=['-h', '--help']))
@click.option('-n', '--name',required=True,
              type=click.Path(),
              help="The annotation file generate by emapper")
@click.option("-f", "--fasta",
              required=True,
              type=click.Path(),
              help="The fasta file contain all of the gene name")
@click.option("-o", "--out",
              required=True,
              type=click.Path(),
              help="The out put dir for the result")

#################################################################################
#主函数（eggonog文件，fa文件，输出文件）
def main(name, fasta, out):
    """
    Parse the eggNOG result and split them to GO KEGG
    """


    #调用函数
    egg_title, egg_dict = parse_eggnog(name)
    
    #BIOPYTHON包读取fasta头部序列id，创建genes_id_list
    genes_id_list = []
    for record in SeqIO.parse(fasta, "fasta"):
        genes_id_list.append(record.id)


    #eggNOG汇总
    dir_eggnog = os.path.join(out, "eggNOG")
    os.makedirs(dir_eggnog, exist_ok=True)
    file_eggnog = os.path.join(dir_eggnog, "eggNOG.xls")
    
    #创建eggNOG.xls，写成excel
    with open(file_eggnog, 'w') as OUT:
        #写入标题
        print(*egg_title, sep="\t", file=OUT)
        #写入内容
        for i, j in egg_dict.items():
            print(*[i, *j], sep="\t", file=OUT)

    #GO
    dir_go = os.path.join(out, "GO")
    os.makedirs(dir_go, exist_ok=True)
    file_go = os.path.join(dir_go, "all.annot")
    
    with open(file_go, 'w') as OUT:
        for i, j in egg_dict.items():
            #第9列COs号，分号间隔
            if j[8] != "-":
                for k in j[8].strip().split(','):
                    print(*[i, k], sep="\t", file=OUT)

    # KEGG
    dir_kegg = os.path.join(out, "KEGG")
    os.makedirs(dir_kegg, exist_ok=True)
    file_kegg = os.path.join(dir_kegg, "all.kopath")
    
    with open(file_kegg, 'w') as OUT:
        print("###version 4.0", file=OUT)
        for gene_id in genes_id_list:
            if gene_id in egg_dict:
                
                #录入第12列ko号。第12列pathway
                KO      =  egg_dict[gene_id][10]   
                pathway =  egg_dict[gene_id][11]
                
                if KO == "-":
                    print(gene_id, file=OUT)
                
                else:
                    #KO格式为"ko:K07733",不需要字符"ko:"，需截取前三个字符串
                    KO = KO[3:]
                    if pathway == '-':
                        print(*[gene_id, KO], sep="\t", file=OUT)
                    
                    else:
                        #pathway格式 "ko00260,ko00270..."
                        #目的截取数值，正则，提取的type为列表，再用join转为字符串
                        pathways = ','.join(set(re.findall("(\d+)", pathway)))
                        print(*[gene_id, KO, pathways], sep="\t", file=OUT)
            
            else:
                print(gene_id, file=OUT)

    # COG分类
    dir_cog = os.path.join(out, "COG")
    os.makedirs(dir_cog, exist_ok=True)
    file_cog = os.path.join(dir_cog, "all.COG.class.xls")

    #COG颜色配置"info_cog"
    #info_cog 格式 "{'A': "#FF9900\tRNA processing and modification",.....}"
    
    res_cog = {i: 0 for i in info_cog.keys()}   ##创建cog分类字典{'A': 0, 'B': 0, 'C': 0...}
    
    for i, j in egg_dict.items():
        if j[5] == '-':
            pass
        
        else:
            #计数
            for tag in j[5]:
                res_cog[tag] += 1

    with open(file_cog, 'w') as OUT:
        print(*["Code", "FunctionalCategories", "GeneNumber", "ColorCode"], sep="\t", file=OUT)
        
        for code in info_cog.keys():
            color, cate = info_cog[code].strip().split("\t")
            arr = [code, cate, res_cog[code], color]
            print(*arr, sep="\t", file=OUT)


if __name__ == "__main__":
    main()