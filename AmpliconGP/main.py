#!/usr/bin/env python

# @CreateTime       : 2022/05/19
# @Author           : mengxf
# @version          : v1.1
# @LastModified     : 2022/12/29

import argparse
import logging
from lib import AGP


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

def get_args():
    parser = argparse.ArgumentParser(
        description="二代测序通用扩增子分析流程 - 多样本并行分析."
    )
    parser.add_argument("-i", "--inyaml", required=True, help="输入信息YAML文件.")
    parser.add_argument("-d", "--dryrun", action="store_true", help="是否运行程序,如选择生成shell脚本不运行.")
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()
    pipe = AGP.Analysis(inyaml=args.inyaml)
    pipe.execute(dryrun=args.dryrun)
    