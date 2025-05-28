from asyncore import read
from email import header
import os
import sys
import pandas as pd

gene_mark_file=sys.argv[1]
outfile=sys.argv[2]


'''
GFF文件是以tab键分割的9列组成，以下为每一列的对应信息：
seq_id：序列的编号，一般为chr或者scanfold编号；
source: 注释的来源，一般为数据库或者注释的机构，如果未知，则用点“.”代替
type: 注释信息的类型，比如Gene、cDNA、mRNA、CDS等
start: 该基因或转录本在参考序列上的起始位置；(从1开始，包含);
end: 该基因或转录本在参考序列上的终止位置；(从1开始，包含);
score: 得分，数字，是注释信息可能性的说明，可以是序列相似性比对时的E-values值或者基因预测是的P-values值，.表示为空;
strand: 该基因或转录本位于参考序列的正链(+)或负链(-)上;
phase: 仅对注释类型为“CDS”有效，表示起始编码的位置，有效值为0、12. (对于编码蛋白质的CDS来说，本列指定下一个密码子开始的位置。每3个核苷酸翻译一个氨基酸，从0开始，CDS的起始位置，除以3，余数就是这个值，，表示到达下一个密码子需要跳过的碱基个数。该编码区第一个密码子的位置，取值0,1,2。0表示该编码框的第一个密码子第一个碱基位于其5’末端；1表示该编码框的第一个密码子的第一个碱基位于该编码区外；2表示该编码框的第一个密码子的第一、二个碱基位于该编码区外；如果Feature为CDS时，必须指明具体值。)；

'''

with open(gene_mark_file,'r')as f,open(outfile,'w')as OUT:
    di={}
    print("基因类型\t数量",file=OUT)

    for line in f:
        if line.startswith("#"):
            continue
        else:
            gene_type=line.strip().split('\t')[2]
            
            if gene_type not in di.keys():
                di[gene_type]=int(1)
                
            else:
                di[gene_type]+=1
    
    for key in di.keys():
        print(f"{key}\t{di[key]}",file=OUT)


