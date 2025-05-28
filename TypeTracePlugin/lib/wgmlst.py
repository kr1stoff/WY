#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2024/2/28 14:25
# @Last Modified by:   Ming
# @Last Modified time: 2024/3/5 14:36
import logging
from pathlib import Path
from subprocess import run

from .cgmlst import cgMLST

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class wgMLST(cgMLST):
    """
    wgMLST 流程
    """

    def run(self):
        """
        生成脚本并执行
        """
        blast_index = Path(f"{self.database['wgmlst']}/{self.subject}/{self.subject}")
        f_scheme = Path(f"{self.database['wgmlst']}/{self.subject}/{self.subject}.txt")
        f_length = Path(f"{self.database['wgmlst']}/{self.subject}/{self.subject}.length")

        self.command = f"""# wgMLST
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
        cmd = f"""cp -rf {self.out}/cgMLST.match.info.tsv {d_upload}/wgMLST.match.info.tsv
if [ -f {self.out}/cgMLST.stat.tsv ]; then
  sed s,cgST,wgST,g -i {self.out}/cgMLST.stat.tsv
  cp -rf {self.out}/cgMLST.stat.tsv {d_upload}/wgMLST.stat.tsv
fi"""
        run(cmd, shell=True, check=True)

    @property
    def result(self):
        """
        获取wgmlst分型结果
        """
        f_result = Path(f"{self.out}/cgMLST.stat.tsv").absolute()
        if f_result.exists():
            with open(f_result, 'r') as IN:
                for line in IN:
                    arr = line.strip().split("\t")
                    if arr[0] == "wgST" or arr[0] == "最相似wgST":
                        return arr[1]
            return None
        else:
            return None
