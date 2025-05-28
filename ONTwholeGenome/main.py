#!/usr/bin/env python

# @CreateTime       : 2023/11/29
# @Author           : mengxf
# @version          : v1.4
# @LastModified     : 2024/3/11
# @description      : ONT平台全基因组测序分析流程

import argparse
import logging
from lib.pipe import Pipe, PipeV, PipeVM
from lib.fastq import ONTFastq


# 运行日志
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

# 命令行参数
parser = argparse.ArgumentParser(description=f"ONT平台全基因组测序分析流程.")
parser.add_argument("-i", "--inyaml", required=True, help="输入YAML文件, SOLAR生成.")
parser.add_argument("-o", "--analysis", default=".",
                    help="结果分析目录, 在该目录下根据task_id新建目录保存结果. (default: .)")
parser.add_argument("-w", "--biwf", default="batassy", choices=["batassy", "virassy", "virmap"],
                    help=f"生信工作流(Bioinfomatics Workflow)名称. "
                    f"batassy:细菌组装,virassy:病毒组装,virmap:病毒有参. (default: batassy)")
parser.add_argument("-c", "--conda_activate", default="/home/earthtest/miniconda3/bin/activate",
                    help="Conda activate 软件路径, conda env 中需要包含流程所需的环境. "
                         "(default: /home/earthtest/miniconda3/bin/activate)")

args = parser.parse_args()
inyaml = args.inyaml
analysis = args.analysis
biwf = args.biwf
cond_act = args.conda_activate

# 主流程
if biwf == "virassy":
    # 病毒组装
    pipe = PipeV(inyaml, analysis, cond_act)
elif biwf == "virmap":
    # 病毒拼接
    pipe = PipeVM(inyaml, analysis, cond_act)
else:
    # 细菌组装
    pipe = Pipe(inyaml, analysis, cond_act)

ontfq = ONTFastq(
    indict=pipe.dict_inyaml,
    workdir=pipe.dir_work)

pipe.preparation()
ontfq.copy_fastq()
pipe.nextflow()
pipe.upload()
pipe.report()
pipe.myzip()
