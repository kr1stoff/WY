#!/usr/bin/env python

import os
import logging
import zipfile


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def link_exist(src, dst):
    """
    只复制存在的文件, 并且没创建过
    [src] sourse源文件
    [dst] destination目的文件
    """
    if os.path.isfile(src) and not os.path.exists(dst):
        os.symlink(src, dst)
    else:
        logging.debug(f"不存在文件<{src}> 或 链接已存在<{dst}>")

class MyLoggingInfo():
    """让logging.info带上流程顺序自增标号"""
    def __init__(self):
        self.count = 0
    
    def info(self, loginfo):
        self.count += 1
        logging.info(f"{self.count}.{loginfo}")


def zip_dir(src, dst):
    """
    zip压缩结果文件夹
    参数
        src     要压缩的文件夹
        dst     压缩后的文件
    """
    zip = zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(src):
        arcpath = root.replace(os.path.dirname(src), "")
        for file in files:
            zip.write(os.path.join(root, file), arcname=os.path.join(arcpath, file))
    zip.close()
