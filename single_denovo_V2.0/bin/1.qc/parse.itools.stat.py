#!/usr/bin/python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/7/26

import logging
import os,re
import click
import sys

#运行日志格式设置：时间+运行提示日志
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
help = dict(help_option_names=['-h', '--help'])



#参数
@click.command(context_settings=help)
@click.option('-i', '--infile',required=True,type=click.Path(),help="iTools Fqtools stat result")
@click.option('-m','--mode',required=True,type=click.Choice(['PE','SE']),default='PE',show_default=True,help="测序模式选择,单/双端")




def main(infile,mode):
    logging.info(f"Parse itools stat {infile}")
    infile = os.path.abspath(infile)
    dir = os.path.dirname(infile)

    with open(infile,"r",encoding='utf-8') as f:
        ReadNum,BaseNum,ReadLeng=[],[],[]
        Nbase,gc=[],[]
        q20,q30=[],[]

        for line in f:
            if not line: continue   #跳过空行
            if not line.strip().startswith("#"): continue   #跳过非"#"开头的行

            #获取总reads数/总碱基数/序列长度
            if line.strip().startswith("#ReadNum"): 
                tmp=re.split("[:\t\n]",line)
                ReadNum.append(tmp[1].replace(' ',''))
                BaseNum.append(tmp[3].replace(' ',''))
                ReadLeng.append(tmp[5].replace(' ',''))

            #获取GC%
            if line.strip().startswith("#GC%"):
                tmp=re.split("[:\t\n]",line)
                gc.append(tmp[1].replace(' ',''))
            
            #获取N碱基
            if line.strip().startswith("#N"):
                tmp=re.split("[:\t\n]",line)
                Nbase.append(tmp[2].replace(' ',''))

            #获取Q20/Q30
            if line.strip().startswith("#BaseQ:10--20"):
                tmp=re.split("[:\t\n]",line)
                q20.append(tmp[4].replace(' ',''))

            if line.strip().startswith("#BaseQ:20--30"):
                tmp=re.split("[:\t\n]",line)
                q30.append(tmp[4].replace(' ',''))

    #写结果
    logging.info(f"Output {dir}/cutfq.stat.txt")
    OUT=open(f"{dir}/cutfq.stat.txt","w",encoding='utf-8')
    OUT.write(f"序列\t总reads数\t总碱基数\t每条reads长度\tGC含量\tN碱基含量\tq20\tq30\n")
    OUT.write(f"抽取read1\t{format(int(ReadNum[0]),',')}\t{format(int(BaseNum[0]),',')}\t{ReadLeng[0]}\t{gc[0]}\t{Nbase[0]}\t{q20[0]}\t{q30[0]}\n")

    if mode == 'PE':
        OUT.write(f"抽取read2\t{format(int(ReadNum[1]),',')}\t{format(int(BaseNum[1]),',')}\t{ReadLeng[1]}\t{gc[1]}\t{Nbase[1]}\t{q20[1]}\t{q30[1]}\n")
    OUT.close()


#调用全局函数
if __name__ == "__main__":
    main()