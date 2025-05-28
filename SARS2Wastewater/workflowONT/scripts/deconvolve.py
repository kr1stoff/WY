import re
from collections import defaultdict
import pandas as pd


def parse_header_ANN(annline):
    """return ANN position dict. Allele, Annotation, Annotation_Impact, Gene_Name..."""
    pattern = re.compile("Functional annotations: \'(.*?)\' \">", re.S)
    annres = re.findall(pattern, annline)[0]
    hdrs = annres.strip().split(" | ")
    dicannpos = {en[1]: en[0] for en in enumerate(hdrs)}
    return dicannpos

def detail_line(line, dichd):
    """vcf和snpeff_vcf重复的步骤"""
    linelist = line.strip().split("\t")
    chrom = linelist[dichd["#CHROM"]]
    pos = linelist[dichd["POS"]]
    ref = linelist[dichd["REF"]]
    alt = linelist[dichd["ALT"]]
    info = linelist[dichd["INFO"]]
    dicinfo = {}
    infs = info.strip().split(";")
    for il in infs:
        ils = il.strip().split("=")
        dicinfo[ils[0]] = ils[1]
    dp = dicinfo['DP']
    af = format(float(dicinfo['AF']), '.2%')
    adp = str(sum([int(i) for i in dicinfo['DP4'].split(',')[-2:]])) #DP4=1,25,550,684
    return chrom, pos, ref, alt, adp, dp, af, dicinfo

def get_mutation_label(ANN, dicannpos, pos, ref, alt):
    """
    在INFO列提取HGVS注释信息 ''.
    返回 [Gene_Name, HGVS.c, HGVS.p]
    220927 DNA变异和氨基酸变异要按照 nextclade 写法
    """
    anns = ANN.strip().split(",")
    annitems = anns[0].strip().split("|")
    gnm = annitems[dicannpos["Gene_Name"]]
    hgvsp = annitems[dicannpos["HGVS.p"]]
    pvar = hgvsp.replace("p.", "")
    pvar = change_3to1_letter(pvar)
    #核酸变异, 按照usher_barcodes的写法
    nvar = f'{ref}{pos}{alt}'
    #是否注释到氨基酸变异
    if pvar:
        annomut = f'{gnm}:{pvar}'
    else:
        if len(ref) == len(alt):
            annomut = f'NUC:{ref}{pos}{alt}'
        elif len(ref) > len(alt):
            annomut = f'DEL:{pos}:{len(ref)-len(alt)}'
        elif len(ref) < len(alt):
            annomut = f'INS:{pos}:{len(alt)-len(ref)}'
    return nvar, annomut

def change_3to1_letter(pvar):
    """获取氨基酸缩写三字母和单字母对照, 将三字母的hgvsp转成单字母nextclade变异形式"""
    dicaa = {
        'Gly': 'G',
        'Ala': 'A',
        'Val': 'V',
        'Leu': 'L',
        'Ile': 'I',
        'Pro': 'P',
        'Phe': 'F',
        'Tyr': 'Y',
        'Trp': 'W',
        'Ser': 'S',
        'Thr': 'T',
        'Cys': 'C',
        'Met': 'M',
        'Asn': 'N',
        'Gln': 'Q',
        'Asp': 'D',
        'Glu': 'E',
        'Lys': 'K',
        'Arg': 'R',
        'His': 'H'
    }
    for ltr3 in dicaa:
        if ltr3 in pvar:
            pvar = pvar.replace(ltr3, dicaa[ltr3])
    return pvar

def get_top_lineages_mutations(topnb=5):
    """
    获取高丰度分型及分型对应的变异, 生成分型-突变字典

    defaultdict(list,
                {'A11201G': ['AY.100'],
                'A23403G': ['AY.100', 'XBB.2.6', 'BF.7.15', 'Q.1', 'Q.2'],
                'C10029T': ['AY.100', 'XBB.2.6', 'BF.7.15']})
    """
    dictoplngsmuts = defaultdict(list)
    #高丰度分型
    dfla = pd.read_table(lngsabdc, sep='\t')
    toplngs = dfla.iloc[:topnb,]['分型'].to_list()
    #freyja user_barcodes 突变-分型对应表
    dfub = pd.read_csv(usherbrcd, sep=',', index_col=0)
    #生成分型-突变字典
    for lng in toplngs:
        #不在user_barcode里面的分型,跳过. 理论上没有,usher_barcode和丰度结果都来自freyja,保持版本一致  
        if lng not in dfub.index:
            continue
        dflm = dfub.loc[lng,].T
        muts = dflm[dflm != 0.0].index.to_list()
        for mut in muts:
            dictoplngsmuts[mut].append(lng)
    return dictoplngsmuts

def vcf2table() -> None:
    with open(invcf) as f,  open(deconvolvetbl, 'w') as g, open(displaytbl, 'w') as h:
        g.write('参考基因组\t位置\t参考碱基\t替换碱基\t突变频率\t突变深度\t总深度\t突变名称\t混合谱系\n')
        h.write('参考基因组\t位置\t参考碱基\t替换碱基\t突变频率\t突变深度\t总深度\t突变名称\t混合谱系\n')
        for line in f:
            if line.startswith("##"):
                if "ID=ANN" in line:
                    dicannpos = parse_header_ANN(line)
            elif line.startswith("#CHROM"):
                hds = line.strip().split("\t")
                dichd = {hie[1]: hie[0] for hie in enumerate(hds)}
            else:
                chrom, pos, ref, alt, adp, dp, af, dicinfo = detail_line(line, dichd)
                nvar, annomut = get_mutation_label(dicinfo["ANN"], dicannpos, pos, ref, alt)
                dictoplngsmuts = get_top_lineages_mutations()
                if nvar in dictoplngsmuts:
                    lineages = ','.join(dictoplngsmuts[nvar])
                else:
                    lineages = 'None found'
                #写入输出表格.deconvolvetbl全展示, displaytbl仅展示不包含"None found"
                g.write('\t'.join([chrom, pos, ref, alt, af, adp, dp, annomut, lineages]) + '\n')
                if lineages != "None found":
                    h.write('\t'.join([chrom, pos, ref, alt, af, adp, dp, annomut, lineages]) + '\n')

if __name__ == '__main__':
    invcf = snakemake.input[0]
    lngsabdc = snakemake.input[1]
    usherbrcd = snakemake.params[0]
    deconvolvetbl = snakemake.output[0]
    displaytbl = snakemake.output[1]
    vcf2table()
