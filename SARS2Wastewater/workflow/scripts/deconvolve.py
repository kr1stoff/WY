from collections import defaultdict
import csv
import pandas as pd
import logging
import pdb

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

def deconvolve_display() -> None:
    with open(ivartbl, encoding='cp1252') as f, open(deconvolvetbl, 'w') as g, open(displaytbl, 'w') as h:
        g.write('参考基因组\t位置\t参考碱基\t替换碱基\t突变频率\t突变深度\t总深度\t突变名称\t混合谱系\n')
        h.write('参考基因组\t位置\t参考碱基\t替换碱基\t突变频率\t突变深度\t总深度\t突变名称\t混合谱系\n')
        reader = csv.reader(f, delimiter='\t')
        #表头字典
        header = f.readline().strip().split('\t')
        dichdidx = {itm:idx for idx,itm in enumerate(header)}
        #迭代行
        prevmutpos = 0
        for row in reader:
            PVAL = float(row[dichdidx['PVAL']])
            POS = int(row[dichdidx['POS']])
            # Parse mutation data from file
            REF = row[dichdidx['REF']]
            ALT = row[dichdidx['ALT']]
            ALT_FREQ = float(row[dichdidx['ALT_FREQ']])
            GFF_FEATURE = row[dichdidx['GFF_FEATURE']]
            REF_AA = row[dichdidx['REF_AA']]
            ALT_AA = row[dichdidx['ALT_AA']]
            POS_AA = row[dichdidx['POS_AA']]
            REGION = row[dichdidx['REGION']]
            ALT_DP = row[dichdidx['ALT_DP']]
            TOTAL_DP = row[dichdidx['TOTAL_DP']]
            #跳过: 1.p-value阈值0.01; 2.突变位置>29903; 3.突变频率小于5%; 4.重复位置突变跳过
            if (PVAL > 0.01) or (POS > 29903) or (ALT_FREQ < 0.05) or (prevmutpos == POS):
                continue
            #刷新"上一个"突变位置
            prevmutpos = POS
            #突变注释 NUC:T12345C, DEL:22222:9, ORF1AB:K47R
            if GFF_FEATURE != 'NA': #注释到蛋白
                annomut = f'{GFF_FEATURE.split(":")[0]}:{REF_AA}{POS_AA}{REF_AA}'
            elif (GFF_FEATURE == 'NA') and (len(REF) != len(ALT)): #没注释到蛋白的DEL
                annomut = f'DEL:{POS}:{len(ALT)-1}' #去掉"-"负号的长度
            elif (GFF_FEATURE == 'NA') and (len(REF) == len(ALT)): #没注释到蛋白的SNP
                annomut = f'NUC:{REF}{POS}{ALT}'
            else:
                logging.warning('未定义的突变类型.')
                continue
            #突变对应分型
            nucmut = f'{REF}{POS}{ALT}'
            if nucmut in dictoplngsmuts:
                lineages = ','.join(dictoplngsmuts[nucmut])
            else:
                lineages = 'None found'
            #写入输出表格.deconvolvetbl全展示, displaytbl仅展示不包含"None found"
            g.write('\t'.join([REGION, str(POS), REF, ALT, f'{ALT_FREQ:.2%}', ALT_DP, TOTAL_DP, annomut, lineages]) + '\n')
            if lineages != "None found":
                h.write('\t'.join([REGION, str(POS), REF, ALT, f'{ALT_FREQ:.2%}', ALT_DP, TOTAL_DP, annomut, lineages]) + '\n')

if __name__ == '__main__':
    ivartbl = snakemake.input[0]
    lngsabdc = snakemake.input[1]
    deconvolvetbl = snakemake.output[0]
    displaytbl = snakemake.output[1]
    usherbrcd = snakemake.params[0]
    dictoplngsmuts = get_top_lineages_mutations()
    deconvolve_display()
