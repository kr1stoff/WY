#!/usr/bin/python3
#-*- coding:utf-8 -*-
import os,re,sys,logging
import gzip
import click
from Bio import SeqIO
from Bio.SeqUtils import GC
import sys

logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


def bwa(id,fq1,fq2,out,ref):
    umsi_id = str(os.path.basename(ref)).split('.')[0]
    sam = f"{out}/{id}_{umsi_id}.sam"
    with open(f"{out}/work.sh" ,'a' ,encoding='utf-8') as w:
        if fq2:
            w.write(f"bwa mem -t 12 -Y {ref} {fq1} {fq2} -o {sam} \n")
        else:
            w.write(f"bwa mem -t 12 -Y {ref} {fq1} -o {sam} \n")

    os.system(f"bash {out}/work.sh")
    return sam


def read_sam(sample_id,sam,out):
    total_reads = 0
    umsi_dic = {}
    with open(sam,'r',encoding='utf-8') as f:
        for line in f:
            if line.strip().startswith("@"):continue
            dict=line.strip().split('\t')
            id,flag,ID_spike_id,cigar,seq,NM_info = dict[0],dict[1],dict[2],dict[5],dict[9],dict[11]
            total_reads += 1
            #跳过比对不上的
            if flag==4:continue

            #跳过cigar = *
            if re.search("[^\dMIDS]",cigar):  continue 

            cigar_num=re.split("M|I|D|S",cigar)
            cigar_str=re.split("\d+",cigar)

            #构造格式：43M7S ['43', '7'] ['M', 'S']
            del cigar_num[-1]
            del cigar_str[0]

            for tmp_str,tmp_num in zip(cigar_str,cigar_num):
                #记录该reads比对上的碱基数
                if tmp_str=='M':
                    map_num=int(tmp_num)
            
            #match匹配率
            if map_num/len(seq) <0.9: continue  

            NM=int(NM_info.split(':')[-1])      #NM:i:0

            #mismatch错配率
            if NM/map_num >0.04: continue

            if ID_spike_id not in umsi_dic.keys():
                umsi_dic.setdefault(ID_spike_id,0)
            
            umsi_dic[ID_spike_id] += 1

    
    #追加写入结果
    with open(f"{out}/umsi.stat.txt" ,'a',encoding='utf-8')as w:
        for umsi_id,values in umsi_dic.items():
            w.write(f"{sample_id}\t{umsi_id}\t{total_reads}\t{values}\n")



##########################################################################################
@click.command(context_settings = dict(help_option_names=['-h', '--help']))
@click.option('--out','-o', type=click.Path(),help="Out dir",default='.',show_default=True)
@click.option('-ref', type=click.Path(),help="UMSI fa")
@click.option('-fq1', type=click.Path(),help="fastq1")
@click.option('-fq2', required=False,type=click.Path(),help="fastq2")
@click.option('-id', type=click.STRING,help="sample id")



def main(id ,fq1 ,fq2 ,out ,ref):
    """
    适用于单双端fq，1个样本可含多个内参样本
    """
    sam = bwa(id ,fq1 ,fq2 ,out ,ref)
    read_sam(id,sam,out)
    



if __name__ == '__main__':
    main()