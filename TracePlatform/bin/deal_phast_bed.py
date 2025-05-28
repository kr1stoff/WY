#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/12/15 15:25
# @Last Modified by:   Ming
# @Last Modified time: 2022/12/15 15:25
import logging
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Function
def is_merge(region1, region2, distance):
    """
    是否合并区间

    :param region1: 区间1的信息
    :param region2: 区间2的信息
    :param distance: 区间距离阈值
    """
    if (region2[1] + 1) - (region1[2] - 1) <= distance:
        return True
    else:
        return False


def print_large_region(region, flag, region_min, merge_max, name, f_handle):
    """
    打印大区域

    :param region: The region bed info
    :param flag: The flag num for the fourth column
    :param region_min: The min length of a region
    :param merge_max: The max length of a region
    :param f_handle: The output file handle
    """
    start = region[1]
    end = min(region[1] + merge_max, region[2])
    while end <= region[2] and start != region[2]:
        # 判断切分后剩余的区域是否小于region_min
        if region[2] - end + 1 > region_min:
            print(*[region[0], start, end, f"{name}-phast{flag}"], sep="\t", file=f_handle)
        else:
            end = region[2]
            print(*[region[0], start, end, f"{name}-phast{flag}"], sep="\t", file=f_handle)
        start = end
        end = min(end + merge_max, region[2])
        flag += 1


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--phast',
              required=True,
              type=click.Path(),
              help="The phast dir of the pipe")
@click.option('--name',
              required=True,
              type=click.Path(),
              help="The Latin name")
@click.option('--merge',
              default=50,
              show_default=True,
              type=int,
              help="当区间距离小于此值时合并区间")
@click.option('--region_min',
              default=200,
              show_default=True,
              type=int,
              help="区间的最小值")
@click.option('--merge_max',
              default=2000,
              show_default=True,
              type=int,
              help="所合并区间的最大长度")
@click.option('--out',
              required=True,
              type=click.Path(),
              help="The out put bed file name")
def main(phast, name, merge, region_min, merge_max, out):
    """
    合并并处理PhastCons生成的结果文件
    """
    d_phast = Path(phast).absolute()
    logger.info(f"Start to Deal the Phast result")
    f_out = Path(out).absolute()

    flag = 1
    with open(f_out, 'w') as OUT:
        f_in = d_phast.joinpath(f"mostcons.bed")
        last_region = None
        with open(f_in, 'r') as IN:
            for line in IN:
                arr = line.strip().split("\t")
                arr[1] = int(arr[1])
                arr[2] = int(arr[2])
                if merge:
                    ## 合并
                    if last_region is None:
                        last_region = arr
                        continue
                    else:
                        if is_merge(last_region, arr, merge):
                            last_region[2] = arr[2]
                        else:
                            if last_region[2] - last_region[1] > merge_max:
                                print_large_region(last_region, flag, region_min, merge_max, name, OUT)
                            else:
                                print(*[last_region[0], last_region[1], last_region[2], f"{name}-phast{flag}"],
                                      sep="\t", file=OUT)
                            last_region = arr
                            flag += 1
                        continue
                else:
                    # 不合并
                    if arr[2] - arr[1] < region_min:
                        print(*[arr[0], arr[1], arr[2], f"{name}-phast{flag}"], sep="\t", file=OUT)
                        flag += 1
                    else:
                        pass
        if last_region:
            if last_region[2] - last_region[1] > merge_max:
                print_large_region(last_region, flag, region_min, merge_max, name, OUT)
            else:
                print(*[last_region[0], last_region[1], last_region[2], f"{name}-phast{flag}"],
                      sep="\t", file=OUT)


if __name__ == "__main__":
    main()
