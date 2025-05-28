#!/usr/bin/env python

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pdb


def get_change_len(genome_len=29903, similarity=0.4):
    """
    根据基因组大小和相似读,算出需要改变的基因组长度
    : [genome_len]  基因组长度
    : [similarity]  相似度
    """
    change_len = round(genome_len - similarity * genome_len)
    return change_len

def simu(fa, outfile, genome_len, similarity):
    record = list(SeqIO.parse(fa, "fasta").records)[0]
    seq = record.seq
    ATGC = list("ATGC")
    change_len_40 = get_change_len(genome_len=genome_len, similarity=similarity)
    list_index = list(range(29903))
    list_variants_index = random.sample(list_index, k=change_len_40)    # 随机抽出发生突变的位置
    for vid in list_variants_index:
        ref = seq[vid]
        ATGC_copy = ATGC.copy()    # 先拷贝一下,remove原地改变列表
        ATGC_copy.remove(ref)   # 外面操作remove
        alt = random.sample(ATGC_copy, k=1)[0]
        seq = seq[:vid] + alt + seq[vid+1:]
    out_records = [SeqRecord(seq, id=f"sars2_{int(similarity*100)}", description="")]
    SeqIO.write(out_records, outfile, "fasta")

# 输入,获取新冠序列
fa = "/sdbb/bioinfor/Database/references/NC_045512.2/NC_045512.2.fasta"
# 输出
for simi in [40, 60, 80]:
    outfile = f"/sdbb/bioinfor/mengxf/Project/6.simulate/results/220511-sars2/sars2_{simi}.fa"
    similarity = simi / 100
    simu(fa, outfile, 29903, similarity)

"""# 输出
dict_variants = {"snp": [], "del": [], "ins": []} # 生成了哪些变异记录下来

# 输入,获取新冠序列
fa = "/sdbb/bioinfor/Database/references/NC_045512.2/NC_045512.2.fasta"
record = list(SeqIO.parse(fa, "fasta").records)[0]
seq = record.seq

# 上面那株奥密克戎影响78bp, indels影响30bp, del与ins 4:1
change_len_40 = get_change_len(genome_len=29903, similarity=0.4)
indels_len_40 = round(change_len_40 * 30 / 78)
del_len_40 = round(indels_len_40 * 4 / 5)
ins_len_40 = round(indels_len_40 * 1 / 5)
snp_len_40 = change_len_40 - del_len_40 - ins_len_40

# 下标列表
now_list_index = list(range(29903)) # 下标列表,抽一次对应去掉部分下标

# del
now_del_len_40 = del_len_40 # 抽一次长度减一次,减到不大于10为止
while now_del_len_40 > 10:
    var_length = random.randint(2, 50) # 暂时模拟indel突变长度在2-50bp
    now_del_len_40 -= var_length
    var_index = random.sample(now_list_index, k=1)[0]
    dict_variants["del"].append(f"{var_index}-{var_index+var_length}")
    for rl in range(var_length):
        print(var_index+rl)
        now_list_index.remove(var_index+rl)

pdb.set_trace()"""
