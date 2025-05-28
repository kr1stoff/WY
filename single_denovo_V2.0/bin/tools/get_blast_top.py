#!/usr/bin/env python
from encodings import utf_8_sig
import os
import sys
import logging
import click
import pandas as pd
import xlwt


logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
help = dict(help_option_names=['-h', '--help'])

#封装程序  main(必须参数，说明)
@click.command(context_settings=help)

@click.option('-i', '--infile',required=True,type=click.Path(),help="input file")
@click.option('-o', '--out',required=True,type=click.Path(),help="output file")
@click.option('-n', '--num',required=True,type=int,help="a number to get top_n")

def main(infile,num,out):
    '''
    According to the first_column
    get the top_n record
    '''
    infile=os.path.abspath(infile)
    out=os.path.abspath(out)

    if os.path.getsize(infile) == int(0) or not os.path.exists(infile):
        open(out,encoding='utf-8').write("比对基因ID\t毒力蛋白ID\t毒力基因简称\t毒力蛋白全称\t微生物体\t覆盖度(%)\t置信度(%)\t参考基因序列长度\t比对基因序列长度\t对齐序列长度\t比对得分\t测序深度")
        dir=os.path.dirname(out)

        open(f"{dir}/top{num}_vir_gene.txt",encoding='utf-8').write("比对基因ID\t毒力蛋白ID\t毒力基因简称\t毒力蛋白全称\t微生物体\t覆盖度(%)\t置信度(%)\t参考基因序列长度\t比对基因序列长度\t对齐序列长度\t比对得分\t测序深度")



    with open(infile,'r') as f,open(out,'w') as OUT:
        logging.info(f"Parse the file {infile}")

        li=set()
        #遍历内容
        for line in f:
            arr = line.strip().split('\t')
            
            #键不在
            if arr[0] not in li:
                print(line.strip(),file=OUT)  

                #初始化变量                 
                flag=int(num)-1
                li.add(arr[0])
            
            #键在
            else:
                if flag != 0:
                    flag-=1
                    print(line.strip(),file=OUT)

    dir=os.path.dirname(out)
    data=pd.read_csv(out,sep="\t",encoding='utf-8')
    pd.DataFrame(data).to_excel(f"{dir}/top{num}_virulence_gene.xlsx",index=False,encoding='utf-8')

    logging.info(f"output file is {out}")   
    

#调用全局函数
if __name__ == "__main__":
    main()
                
        

