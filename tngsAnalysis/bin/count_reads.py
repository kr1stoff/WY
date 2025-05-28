#!/usr/bin/python3
# @auther: zhuzi
# @date: 2022-12-12

import sys
import click
import logging
import re

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)


# 判断：identity是否通过
def isPassIdentity(identity,match,mismatch):
    if (match - mismatch)/match >= identity:
        return True

def filter_sam(sam,identity,out):
    map_reads_Dic = {}
    total_reads = 0
    with open(sam,'r',encoding='utf-8') as f ,open(out,'w',encoding='utf-8') as w:
        for line in f:
            if line.strip().startswith("@"):
                w.write(f"{line}")
                continue
            lines = line.strip().split('\t')
            flag,match_id,cigar = int(lines[1]),lines[2],lines[5]
            if (flag & 256 == 256) or (flag & 2048 == 2048): continue

            total_reads += 1
            if flag & 4 == 4 : continue
            match = sum([int(i.strip('M')) for i in re.findall('\d+M', cigar)])

            mismatch = re.findall('NM:i:(\d+)',line)[0]
            if not mismatch: mismatch = 0
            
            if not isPassIdentity(identity,match,int(mismatch)): continue
            w.write(f"{line}")   # write filter sam

            map_reads_Dic.setdefault(match_id,0)
            map_reads_Dic[match_id] += 1

        return map_reads_Dic,total_reads


# === Option =============================================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('sam_file', type=click.Path())
@click.option('-o','--out',default='./filter.sam',type=click.Path(),help="过滤sam文件")
@click.option('-s','--stat',default='./count_reads.stat',type=click.Path(),help="reads比对统计")
@click.option('--identity',default=0.95,type=click.FLOAT,help="identity，默认[0.95]")

def main(sam_file,out,identity,stat):
    """统计sam文件比对上基因的reads数"""
    map_reads_Dic,total_reads = filter_sam(sam_file,identity,out)
    logging.info(f"[Output]: {out}")
    with open(stat,'w',encoding='utf-8') as w:
        w.write(f"#total_reads:{total_reads}\n")
        w.write(f"#seq\tcount\n")
        for k,v in map_reads_Dic.items():
            w.write(f"{k}\t{v}\n")
    logging.info(f"[Output]: {stat}")


if __name__ == '__main__':
    main()