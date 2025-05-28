#!/usr/bin/env python

import os
import sys
import zipfile


def zip_dir(dir_src, zip_dst):
    """
    压缩目录,和'linux zip -r'一样,不用切目录
        dir_src - 目标文件夹路径
        zip_dst - 目的zip文件路径
    """
    zip = zipfile.ZipFile(zip_dst, "w", zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(dir_src):
        arcpath = root.replace(os.path.dirname(dir_src), "")
        for file in files:
            zip.write(os.path.join(root, file), arcname=os.path.join(arcpath, file))
    zip.close()

if __name__ == "__main__":
    # 读参数
    usage = f"Usage:\n\tpython {os.path.basename(sys.argv[0])} <dir_src> <zip_dst>\n"
    if len(sys.argv) != 3:
        print(usage)
        sys.exit(1)
    dir_src, zip_dst = sys.argv[1:]
    # main
    zip_dir(dir_src, zip_dst)
