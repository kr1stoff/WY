#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Ming
# @Date:   2022-05-29 18:12:37
# @Last Modified by:   Ming
# @Last Modified time: 2022-05-29 18:12:37
import logging
import multiprocessing
import os
from pathlib import Path
from subprocess import run

import yaml

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Some Setting
PATH = os.path.dirname(os.path.abspath(__file__))
CONFIGPATH = os.path.join(PATH, "../config")
BINPATH = os.path.join(PATH, "../bin")
f_software = os.path.join(CONFIGPATH, "software.yml")
f_database = os.path.join(CONFIGPATH, "database.yml")
f_env = os.path.join(CONFIGPATH, "env.yml")


class cgMLST(object):
    """
    cgMLST analysis
    """

    def __init__(self, fasta, subject, out):
        """
        Init the object

        :param fasta: The input genome sequence file
        :param subject: The scheme to use
        :param out: The output directory
        """
        self.env = yaml.safe_load(open(f_env, 'r').read())
        self.software = yaml.safe_load(open(f_software, 'r').read())
        self.database = yaml.safe_load(open(f_database, 'r').read())

        self.subject = f"{self.database['cgmlst']}/{subject}"
        self.subject = subject
        self.fasta = os.path.abspath(fasta)
        self.out = os.path.abspath(out)

        os.makedirs(self.out, exist_ok=True)
        self.script = os.path.join(self.out, "cgmlst.sh")
        self.log = os.path.join(self.out, "cgmlst.log")

        # 机器信息
        self.cpu = multiprocessing.cpu_count()

    def run(self):
        """
        生成脚本并执行
        """
        blast_index = Path(f"{self.database['cgmlst']}/{self.subject}/{self.subject}")
        f_scheme = Path(f"{self.database['cgmlst']}/{self.subject}/{self.subject}.txt")
        f_length = Path(f"{self.database['cgmlst']}/{self.subject}/{self.subject}.length")

        self.command = f"""# cgMLST
set -e
mkdir -p {self.out}/prepare
ln -sf {self.fasta} {self.out}/prepare/Query.fna
if [ -f {f_scheme} ]; then
  {self.software["python"]} {self.software["cgmlst"]}/main.py cgmlst -t 24 --sequence {self.out}/prepare/Query.fna -db {blast_index} --length {f_length} --scheme {f_scheme} -o {self.out}
else
  {self.software["python"]} {self.software["cgmlst"]}/main.py cgmlst -t 24 --sequence {self.out}/prepare/Query.fna -db {blast_index} --length {f_length} -o {self.out}
fi
"""
        with open(self.script, 'w') as OUT:
            print(self.command, file=OUT)
        run(f"{self.software['shell']} {self.script} 2>&1 > {self.log}",
            shell=True, check=True)

    def to_upload(self, d_upload):
        """
        link the result file to dir upload
        """
        cmd = f"""cp -rf {self.out}/cgMLST.match.info.tsv {d_upload}/cgMLST.match.info.tsv
if [ -f {self.out}/cgMLST.stat.tsv ]; then
  cp -rf {self.out}/cgMLST.stat.tsv {d_upload}/
fi"""

        run(cmd, shell=True)

    @property
    def result(self):
        """
        获取cgmlst分型结果
        """
        f_result = Path(f"{self.out}/cgMLST.stat.tsv").absolute()
        if f_result.exists():
            with open(f_result, 'r') as IN:
                for line in IN:
                    arr = line.strip().split("\t")
                    if arr[0] == "cgST" or arr[0] == "最相似cgST":
                        return arr[1]
            return None
        else:
            return None
