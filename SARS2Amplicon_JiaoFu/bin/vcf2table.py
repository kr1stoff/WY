#!/usr/bin/env python
#230112 1.删除变异深度；2.氨基酸变异删掉基因名

import os
import re
import sys
import logging
from pathlib import Path
from openpyxl import Workbook


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

def parse_vcfinfo(info):
    """return INFO dict. DP,AC,AF..."""
    info_dict = dict()
    info_list = info.strip().split(";")
    for il in info_list:
        il_list = il.strip().split("=")
        info_dict[il_list[0]] = il_list[1]
    return info_dict

def parse_header_ANN(annline):
    """return ANN position dict. Allele, Annotation, Annotation_Impact, Gene_Name..."""
    pattern = re.compile("Functional annotations: \'(.*?)\' \">", re.S)
    annres = re.findall(pattern, annline)[0]
    header_list = annres.strip().split(" | ")
    annpos_dict = {en[1]: en[0] for en in enumerate(header_list)}
    return annpos_dict    

def extract_anno_snpeff(ANN, annpos_dict):
    """
    在INFO列提取HGVS注释信息 ''.
    返回 [Gene_Name, HGVS.c, HGVS.p]
    220927 DNA变异和氨基酸变异要按照 nextclade 写法
    """
    ann_list = ANN.strip().split(",")
    # choose position 1 items
    annitem_list = ann_list[0].strip().split("|")
    genename = annitem_list[annpos_dict["Gene_Name"]]
    hgvsc = annitem_list[annpos_dict["HGVS.c"]]
    hgvsp = annitem_list[annpos_dict["HGVS.p"]]
    nvar, pvar = hgvs2nextclade(hgvsc, hgvsp)
    return [genename, nvar, pvar]

def tsv2xlsx(tfile):
    outfile = re.sub(r".txt|.tsv|.xls|.csv", ".xlsx", tfile)
    wb = Workbook()
    ws = wb.active
    with open(tfile, "rt") as fh:
        for line in fh:
            ws.append(line.strip().split("\t"))
    wb.save(outfile)

def get_variant_info(vform, iform):
    """根据VCF最后两列 FORMAR/样本信息 找到测序深度/变异深度/变异频率."""
    vfs = vform.strip().split(":")
    vi = {vi[1]:vi[0] for vi in enumerate(vfs)} # 下标字典
    ifs = iform.strip().split(":")
    ad = list(map(int, ifs[vi["AD"]].strip().split(",")))
    dp = format(int(ifs[vi["DP"]]), ",")
    if len(ad) == 2: # 0/1, 1/1
        ro, _ao = ad # observations of the reference/alternate
        ao = format(_ao, ",")
        af = f"{_ao/int(ifs[vi['DP']])*100:.2f}"
    elif len(ad) == 3: # 1/2 ...
        ro, _ao1, _ao2 = ad
        ao = f"{format(_ao1, ',')};{format(_ao2, ',')}"
        af = f"{_ao1/int(ifs[vi['DP']])*100:.2f};{_ao2/int(ifs[vi['DP']])*100:.2f}"
    else: # 220927 碰到再看看
        logging.warning(ifs[vi["AD"]])
        ao, af = "-", "-"
    return ao, dp, af

def detail_line(line, header_dict):
    """vcf和snpeff_vcf重复的步骤"""
    linelist    = line.strip().split("\t")
    chrom       = linelist[header_dict["#CHROM"]]
    pos         = linelist[header_dict["POS"]]
    ref         = linelist[header_dict["REF"]]
    alt         = linelist[header_dict["ALT"]]
    qual        = float(linelist[header_dict["QUAL"]])
    info        = linelist[header_dict["INFO"]]
    vformat     = linelist[header_dict["FORMAT"]]
    iformat     = linelist[-1]  # FORMAT 对应的信息
    info_dict   = parse_vcfinfo(info)
    ao, dp, af = get_variant_info(vformat, iformat)
    return chrom, pos, ref, alt, qual, ao, dp, af, info, info_dict

def change_3to1_letter(pvar):
    """获取氨基酸缩写三字母和单字母对照, 将三字母的hgvsp转成单字母nextclade变异形式"""
    aacomp = f"{Path(__file__).parents[1]}/etc/amino_acid_comparison.tsv"
    dict_aa = {}
    with open(aacomp, "rt") as fh:
        header = next(fh)
        heads = header.strip().split("\t")
        hi = {it[1]: it[0] for it in enumerate(heads)}
        for line in fh:
            lns = line.strip().split("\t")
            dict_aa[lns[hi["三字母缩写"]]] = lns[hi["单字母符号"]]
    for i in range(2): # 应该只有两个氨基酸吧, 多了就写3
        for tlett in dict_aa: # tlett three letter 三字母字符
            if tlett in pvar:
                pvar = pvar.replace(tlett, dict_aa[tlett])
    return pvar

def hgvs2nextclade(hgvsc, hgvsp):
    """
    hgvs的变异格式改成nextclade的格式 
    c.9601C>T   --> C9601T
    p.Ser135Arg --> S135R
    """
    # 核酸
    nvar = hgvsc.replace("c.", "")
    npat = re.compile("(-?\d+)([A-Z])>([A-Z])")
    if ">" in nvar:
        nres = re.findall(npat, nvar)
        pos, ref, alt = nres[0]
        nvar = "".join([ref, pos, alt])
    # 氨基酸
    pvar = hgvsp.replace("p.", "")
    pvar = change_3to1_letter(pvar)
    return nvar, pvar

def vcf2table(vfile, outfile):
    with open(vfile) as fh, open(outfile, "wt", encoding="utf-8", newline="") as gh:
        gh.write("参考基因组\t变异位置\t参考碱基\t替换碱基\t基因名\tDNA变异\t氨基酸变异\t变异频率%\t测序深度\n")
        for line in fh:
            if line.startswith("##"):
                if "ID=ANN" in line:
                    annpos_dict = parse_header_ANN(line)
                continue
            elif line.startswith("#CHROM"):
                header_list = line.strip().split("\t")
                header_dict = {hie[1]: hie[0] for hie in enumerate(header_list)}
            else:
                chrom, pos, ref, alt, qual, ao, dp, af, info, info_dict = detail_line(line, header_dict)
                ann3_list   = extract_anno_snpeff(info_dict["ANN"], annpos_dict)
                if ("DP" not in info_dict) or (qual < 1): # 1.无DP, 2.质量<1
                    logging.warning(info)
                    continue
                else:
                    outlist = [chrom, pos, ref, alt] + ann3_list + [af, dp]
                    gh.write("\t".join(outlist) + "\n")

def vcf2table_normal(vfile, outfile):
    """相较于vcf2table函数, 处理未注释的VCF文件"""
    with open(vfile) as fh, open(outfile, "wt", encoding="utf-8", newline="") as gh:
        gh.write("参考基因组\t变异位置\t参考碱基\t替换碱基\t变异频率%\t测序深度\n")
        for line in fh:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                header_list = line.strip().split("\t")
                header_dict = {hie[1]: hie[0] for hie in enumerate(header_list)}
            else:
                chrom, pos, ref, alt, qual, ao, dp, af, info, info_dict = detail_line(line, header_dict)
                if ("DP" not in info_dict) or (qual < 1): # 1.无DP, 2.质量<1
                    logging.warning(info)
                    continue
                else:
                    outlist = [chrom, pos, ref, alt] + [af, dp]
                    gh.write("\t".join(outlist) + "\n")

def check_vcf(vfile):
    """查看输入VCF是普通的还是snpEff注释过的"""
    with open(vfile, "rt") as fh:
        context = fh.read()
        if "SnpEff" in context:
            vtype = "snpEff"
        else:
            vtype = "normal"
    return vtype

# main
if __name__ == "__main__":
    if len(sys.argv) != 3 or not os.path.exists(sys.argv[1]):
        sys.exit(f"{Path(sys.argv[0]).name} <in_vcf> <out_tsv>")
    vcf_file, out_tsv = sys.argv[1:]
    vtype = check_vcf(vcf_file)
    if vtype == "snpEff":
        vcf2table(vcf_file, out_tsv)
    elif vtype == "normal":
        vcf2table_normal(vcf_file, out_tsv)
    tsv2xlsx(out_tsv)
    