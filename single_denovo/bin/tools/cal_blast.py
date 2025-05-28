#!/usr/bin/env python3
#-*- coding:utf-8 -*-
import os
import click
import sys
import logging
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


# 区间合并
def merge(record_dic):
    """  {ref_id:{区间：[.....],长度:xxxx}}  """
    new_interval = {}
    for ref_id,intervals_dic in record_dic.items():
        intervals = record_dic[ref_id]['区间']
        intervals.sort(key=lambda x: x[0])  #按sstart进行排序
        merged = []

        # 遍历区间
        for record in intervals:
            # 当前遍历的start位点 > 已经记录的区间
            if not merged or record[0] > merged[-1][-1]: 
                merged.append(record)
            
            # 重叠区处理，哪个大取哪个
            else:               
                merged[-1][-1] = max(merged[-1][-1], record[-1])
        new_interval[ref_id] = merged
    return new_interval



def parse_blast(infile,mismatch_rate,identity_rate):
    
    dir = os.path.dirname(infile)
    record_dic = {}    #可以适配参考基因组不止一条fasta
    filter_mismatch,filter_identity = 0,0 

    with open(infile,'r',encoding='utf-8') as f:
        for lines in f:
            if not lines: continue
            line = lines.strip().split('\t')

            # 错配率筛选mistch_num/align_length
            if int(line[10])/int(line[9]) > float(mismatch_rate):
                filter_mismatch += 1
                continue

            # 置信度筛选idtntity
            if float(line[8]) < float(identity_rate)*100:
                filter_identity += 1
                continue
            
            # 区间处理(sstart send)
            sstart = int(line[6])
            send = int(line[7])
            tmp_interval = [sstart,send] if (sstart < send) else [send,sstart]
            ref_id = line[4]
            if ref_id not in record_dic.keys(): 
                record_dic[ref_id] = {}
                record_dic[ref_id]['区间'] = []
                record_dic[ref_id]['参考基因组长度'] = int(line[5])
                record_dic[ref_id]['对齐长度'] = 0
                record_dic[ref_id]['错配碱基数'] = 0

            record_dic[ref_id]['区间'].append(tmp_interval)

            #置信度相关
            align_length = int(line[9])
            mismatch_num = int(line[10])
            record_dic[ref_id]['对齐长度'] += align_length
            record_dic[ref_id]['错配碱基数'] += mismatch_num



    logging.info(f"Filter {filter_mismatch} record by filter_mismatch")
    logging.info(f"Filter {filter_identity} record by filter_identity")
    logging.info(f"Blast infile 'sseqid' has {len(record_dic.keys())} {record_dic.keys()} record")
    
    # 区间处理
    new_interval = merge(record_dic)

    # 覆盖度/置信度
    OUT = open(f"{dir}/ref_cov_stat.txt",'w')
    OUT.write(f"参考基因组ID\t参考基因组长度\t对齐基因长度\t覆盖度\t置信度\n")
    identity = (record_dic[ref_id]['对齐长度']-record_dic[ref_id]['错配碱基数'])/record_dic[ref_id]['对齐长度']
    for ref_id in new_interval.keys():
        length = 0
        ref_len = record_dic[ref_id]['参考基因组长度']
        for pair in new_interval[ref_id]:
            (start,end) = pair
            length += (end - start + 1)
        cov = length/ref_len
        OUT.write(f"{ref_id}\t{format(ref_len,',')}\t{format(length,',')}\t{cov:.2%}\t{identity:.2%}\n")



##########################################################################################
@click.command(context_settings = dict(help_option_names=['-h', '--help']))
@click.option('--infile','-i', required=True,type=click.Path(),help="blast结果")
@click.option('--mismatch_rate', type=click.FLOAT,default=0.1,show_default=True,help="指定单行结果进行错配率筛选")
@click.option('--identity_rate', type=click.FLOAT,default=0.8,show_default=True,help="指定单行结果进行置信度筛选")


def main(infile,mismatch_rate,identity_rate):
    """
    根据blast结果，进行覆盖度计算，其中可设置一些参数进行balst结果筛选
    blast顺序需满足：-outfmt '6 qseqid qlen qstart qend sseqid slen sstart send pident length mismatch gapopen evalue bitscore'
    """

    infile = os.path.abspath(infile)
    parse_blast(infile,mismatch_rate,identity_rate)



if __name__ == '__main__':
    main()