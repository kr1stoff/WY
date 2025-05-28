#!/usr/bin/env python

import argparse
import logging
import os
import sys
from multiprocessing import Pool
from subprocess import run

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


def execute_linux_commandline(cmd):
    logging.debug("Command Line: " + cmd)
    res = run(cmd, shell=True, encoding="utf-8")
    logging.debug(res.returncode)


def main(sh_file, procs):
    with open(sh_file, "rt", encoding="utf-8") as fh:
        cmd_list = [line.strip() for line in fh]
    # 并行
    pool = Pool(procs)
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
