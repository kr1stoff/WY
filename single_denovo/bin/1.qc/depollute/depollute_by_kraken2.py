#!/usr/bin/env python
# @CreateTime       : 2022/08/22
# @Author           : mengxf
# @LastModified     : 2022/09/22
# 220822 单菌实验原始数据去污染的一种方法,使用kraken2分类,找到优势物种进行污染去除
# 220922 单菌流程允许单端, 加入单端数据处理

import re
import argparse
import logging
import yaml
import shutil
from pathlib import Path
from subprocess import run
from functools import wraps
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys

logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def logit(func):
    """打印运行命令的装饰器"""
    @wraps(func)
    def wrapper(cml):
        logging.debug(f"运行命令: {cml}")
        result = func(cml)
        return result
    return wrapper

@logit
def wraprun(cml):
    """包装run减少重复操作"""
    res = run(cml, shell=True, capture_output=True, encoding="utf-8")
    return res


def get_kraken_optimal(kreport, rank="S"):
    """
    解析kraken2.report 找到最佳比对物种或属等等，默认 种(即筛选报告中rank='S'，中丰度最高的一个)
    返回pandas.Series对象
    """
    ranks = "KPCOFGS"
    if rank not in ranks:
        logging.error("输入分类级别不存在或不是一般分类级别!")
    df = pd.read_table(kreport, sep="\t", header=None, names=["percentage","covered","assigned","rank","taxid","sciname"])
    df["sciname"] = df["sciname"].apply(lambda x: x.strip())
    optimal = df[df["rank"] == rank].sort_values(by=["covered"],ascending=False).iloc[0,]   # 丰度最高的一个taxid
    return optimal 


def get_lineages(taxid):
    """获取所有和目标taxid上下世系的taxids列表"""

    # 获取向下所有子孙taxids列表
    cml = f"{params.taxonkit} list --ids {taxid} --data-dir {params.taxdmp} --indent ''"
    res = wraprun(cml)
    lower_lineages = res.stdout.strip().split("\n")

    # 向上所有祖先taxids
    cml = f"echo {taxid} | {params.taxonkit} lineage -t --data-dir {params.taxdmp}"
    res = wraprun(cml)
    higher_lineages = res.stdout.strip().split("\t")[-1].split(";")[:-1]
    # 整条lineages taxids列表
    return  lower_lineages + higher_lineages

def run_kraken2():
    """运行kraken2 220922 支持单端"""

    # 支持双端
    if len(args.fqs) == 2:
        cml = f"""
{params.kraken2} --threads {params.threads} --paired -db {params.kraken2_db} --output {dir_temp}/{name}.out --report {dir_temp}/{name}.report \
    --classified-out {dir_temp}/{name}.classified#.fq --unclassified-out {dir_temp}/{name}.unclassified#.fq \
    {args.fqs[0]} {args.fqs[1]}
"""

    # 支持单端
    elif len(args.fqs) == 1: 
        cml = f"""
{params.kraken2} --threads {params.threads} -db {params.kraken2_db} --output {dir_temp}/{name}.out --report {dir_temp}/{name}.report \
    --classified-out {dir_temp}/{name}.classified_1.fq --unclassified-out {dir_temp}/{name}.unclassified_1.fq \
    {args.fqs[0]}
"""
    wraprun(cml) # run



def decontaminate():
    """
    使用Biopython 抽去污染后的序列
    220922 支持单端
    """
    def repeat_func1(classified_fq, outfq):
        """重复步骤,抽取双端read"""
        with open(classified_fq, "rt") as ih,\
            open(outfq, "wt", encoding="utf-8", newline="") as gh:
            for title, seq, qual in FastqGeneralIterator(ih):
                try:
                    taxid = re.findall(r"kraken:taxid\|(\d+)", title)[0]
                    if taxid in optimal_lineages:
                        gh.write(f"@{title}\n{seq}\n+\n{qual}\n")
                except Exception as e:
                    logging.warning(e)
                    logging.warning(title)

    repeat_func1(f"{dir_temp}/{name}.classified_1.fq", f"{dir_temp}/{name}.select.1.fq")
    cml = f"""
set -e
cat {dir_temp}/{name}.select.1.fq {dir_temp}/{name}.unclassified_1.fq > {args.prefix}_1.fq
    """
    if len(args.fqs) == 2: # 双端
        repeat_func1(f"{dir_temp}/{name}.classified_2.fq", f"{dir_temp}/{name}.select.2.fq")
        cml += f"""
cat {dir_temp}/{name}.select.2.fq {dir_temp}/{name}.unclassified_2.fq > {args.prefix}_2.fq
        """
    wraprun(cml)

def configuare(yaml_conf):
    """从YAML配置文件获取参数带index数组"""
    dict_params = yaml.safe_load(open(yaml_conf, "rt"))
    params = pd.Series(data=dict_params, index=dict_params.keys())
    return params

def get_args():
    parser = argparse.ArgumentParser()
    parser.description = "去污染程序, 输入单/双端FASTQ, 输出去污染后单/双端FASTQ."
    parser.add_argument("fqs", nargs="+", help="输入单/双FASTQ数据.")
    parser.add_argument("-p", "--prefix", default="./out", help="输出FASTQ前缀,目录加名称. [默认:./out]")
    parser.add_argument("--config", default=Path(__file__).parent.joinpath("decontamination_by_kraken2.yaml"), 
                        help="配置文件.包含软件路径,数据库路径及其他参数.")
    parser.add_argument("--rank", default="S", choices=list("KPCOFGS"), help="过滤使用的分类级别. [默认:S(物种级)]")
    args = parser.parse_args()
    for fil in args.fqs + [args.config]:
        if not Path(fil).exists(): parser.print_help() # 检查不存在,打印帮助
    return args

if __name__ == "__main__":
    args = get_args()
    params = configuare(args.config)
    name = Path(args.prefix).name
    dir_temp = Path(args.prefix).parent.joinpath(f"temp_{name}")
    if Path(dir_temp).exists():
        shutil.rmtree(dir_temp)
    Path(dir_temp).mkdir(exist_ok=True)
    
    # run
    run_kraken2()
    optimal = get_kraken_optimal(f"{dir_temp}/{name}.report", rank=args.rank)
    optimal_lineages = get_lineages(optimal.taxid)
    decontaminate()
    with open(f"{dir_temp}/target.taxid.list",'w',encoding='utf-8') as f:
        print(*optimal_lineages,sep='\n',file=f)
    
