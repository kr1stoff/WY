#!/usr/bin/env python

# @CreateTime       : 2022/05/19
# @Author           : mengxf
# @version          : v1.0
# @LastModified     : 2022/05/26

import logging
import argparse
from lib import ILMN
from lib import MGI


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

def get_args():
    parser = argparse.ArgumentParser(
        description="二代测序Barcode拆分流程, 支持Illumina,MGI, 支持SE/PE, 支持 Single/Dual Index."
    )
    parser.add_argument("-i", "--inyaml", required=True, help="输入信息YAML文件.")
    parser.add_argument("-p", "--platform", required=True, choices=["Illumina", "MGI"], help="测序平台.")
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()
    if args.platform == "Illumina":
        pipe = ILMN.ILMN(inyaml=args.inyaml)
    elif args.platform == "MGI":
        pipe = MGI.MGI(inyaml=args.inyaml)
    pipe.execute()
