#!/usr/bin/env python

# @CreateTime       : 2022/05/16
# @Author           : mengxf
# @version          : v1.1
# @LastModified     : 2022/05/16

# @用意：替代GNU parallel主要功能

import os, sys
import logging
import argparse
from subprocess import run
from multiprocessing import Pool


logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")

def execute_linux_commandline(cmd):
    logging.debug("Command Line: " + cmd)
    res = run(cmd, shell=True, encoding="utf-8")
    logging.debug(res.returncode)

def main(shell_scripts, processes):
    with open(shell_scripts, "rt", encoding="utf-8") as fh:
        cmd_list = [line.strip() for line in fh]
    # 并行
    pool = Pool(processes)
    pool.map_async(execute_linux_commandline, cmd_list)
    logging.info("Waiting for all subprocesses done...")
    pool.close()
    pool.join()
    logging.info("All subprocesses done.")

def get_argparses():
    parser = argparse.ArgumentParser()
    parser.description = "Parallel Program for execute Linux CommandLine!!"
    parser.add_argument("shell_scripts", type=str, help="所有并行脚本的总 shell 文件")
    parser.add_argument("-p", "--processes", default=4, type=int, help="并行的程序数，也可以当作核心数 (default: 4)")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_argparses()
    shell_scripts = args.shell_scripts
    processes = args.processes
    main(shell_scripts, processes)
