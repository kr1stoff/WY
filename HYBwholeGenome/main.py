#!/usr/bin/env python

# @CreateTime       : 2023/12/19
# @Author           : mengxf
# @version          : v1.1
# @LastModified     : 2024/03/11
# @description      : ONT+ILMN混装全基因组测序分析流程

import argparse
import logging
# custom
from lib.pipe import Pipe, PipeV
from lib.fastq import HYBfastq


# 运行日志
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


# 命令行参数
parser = argparse.ArgumentParser(description=f"ONT+ILMN混装全基因组测序分析流程.")
parser.add_argument("-i", "--inyaml", required=True, help="输入YAML文件, SOLAR生成.")
parser.add_argument("-o", "--analysis", default=".",
                    help="结果分析目录, 在该目录下根据task_id新建目录保存结果. (default: .)")
parser.add_argument("-w", "--biwf", default="batassy", choices=["batassy", "virassy"],
                    help=f"生信工作流(Bioinfomatics Workflow)名称. batassy:细菌组装,virassy:病毒组装. (default: batassy)")
parser.add_argument("-c", "--conda_activate", default="/home/earthtest/miniconda3/bin/activate",
                    help="Conda activate 软件路径, conda env 中需要包含流程所需的环境."
                    "(default: /home/earthtest/miniconda3/bin/activate)")

args = parser.parse_args()
inyaml = args.inyaml
analysis = args.analysis
biwf = args.biwf
cond_acti = args.conda_activate

# 主流程
if biwf == "virassy":
    # 病毒组装
    pipe = PipeV(inyaml, analysis, cond_acti)
else:
    # 细菌组装
    pipe = Pipe(inyaml, analysis, cond_acti)

hybfq = HYBfastq(
    indict=pipe.dict_inyml,
    workdir=pipe.dir_work)

pipe.preparation()
hybfq.copy_fastq()
hybfq.generate_samplesheet()
pipe.nextflow()
pipe.upload()
pipe.report()
pipe.myzip()
