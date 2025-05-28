#!/usr/bin/env python
# @CreateTime       : 2023/01/09
# @Author           : mengxf
# @version          : v1.0.3
# @LastModified     : 2023/01/17
# @description      : [交付版本]二代测序新冠扩增子生信分析流程

import argparse
import logging
from lib import SARS2


# 设置运行日志
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

def get_args():
    parser = argparse.ArgumentParser(description="二代测序新冠扩增子分析流程 - 多样本并行分析.")
    parser.add_argument("-n", "--dryrun", action="store_true", help="是否运行程序,如选择生成shell脚本不运行.")
    parser.add_argument("-i", "--inyaml", required=True, help="输入信息YAML文件.")
    parser.add_argument('--trim_software', default='ivar', choices=["ivar","fgbio","ktrim"],
                        help='去引物程序选择. ["ivar","fgbio","ktrim"]')
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()
    pipe = SARS2.SARS2(args)
    pipe.execute(dryrun=args.dryrun)
