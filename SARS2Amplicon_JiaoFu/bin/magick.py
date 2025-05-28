#!/usr/bin/env python
import os
import sys
from glob import glob
from subprocess import run


def magick_dir(dir_fig, path_magick, suffix_from, suffix_to):
    """
    整个目录下的图片转换格式
    [dir_fig]       图片目录
    [path_magick]   magick软件路径
    [suffix_from]   模板图片后缀
    [suffix_to]     需要生成图片后缀
    """
    svgs = glob(f"{dir_fig}/*.{suffix_from}")
    if suffix_from == "svg" and suffix_to == "png":
        for svg in svgs:
            run(f"{path_magick} -density 300 -colorspace RGB {svg} {svg.replace(suffix_from, suffix_to)}", 
            shell=True, check=True)
    else:
        print("mxf: 添加判断.")

if __name__ == "__main__":
    # 获取参数
    usage = f"Usage:\n\tpython {os.path.basename(sys.argv[0])} <magick_path> <figure_dir>\n"
    if len(sys.argv) != 3:
        print(usage)
        sys.exit(1)
    path_magick, dir_figure = sys.argv[1:]
    # main
    magick_dir(dir_figure, path_magick, "svg", "png")
