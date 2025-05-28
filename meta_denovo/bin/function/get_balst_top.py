#!/usr/bin/env python
import os
import sys
import logging
import click

#运行日志格式设置：时间+运行提示日志
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s')

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

    logging.info(f"output file is {out}")   
    

#调用全局函数
if __name__ == "__main__":
    main()
                
        

