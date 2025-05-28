# @CreateTime       : 2022/04/26
# @Author           : mengxf
# @version          : v2.1.1
# @LastModified     : 2023/10/31
# @Describtion      : 溯源进化树分析软件

import sys
import argparse
import logging
from lib.TraceSNP import PhyloSNP
from lib.TraceWGS import PhyloWGS
from lib.TraceCORE import PhyloCORE
from lib import general


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format='%(levelname)s - %(asctime)s - %(message)s',
    datefmt='%Y/%m/%d %H:%M:%S'
)

# 命令行参数
args = general.get_args()

# 判断并运行流程
if args.pipe == "wgs":
    pipe = PhyloWGS(args.inyaml)
elif args.pipe == "snp":
    pipe = PhyloSNP(args.inyaml)
elif args.pipe == "core":
    pipe = PhyloCORE(args.inyaml)
pipe.execute(args.dryrun)
