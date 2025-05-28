#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Ming
# @Date:   2022-05-29 18:12:37
# @Last Modified by:   Ming
# @Last Modified time: 2022-05-29 18:12:37
import logging
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


class MLST(object):
    """
    MLST analysis
    """

    def __init__(self, fasta, subject, out):
        """
        Init the object

        :param fasta: The input fasta sequence
        :param subject: The scheme to use
        :param out: The output dir
        """
        self.env = yaml.safe_load(open(f_env, 'r').read())
        self.software = yaml.safe_load(open(f_software, 'r').read())
        self.database = yaml.safe_load(open(f_database, 'r').read())

        self.subject = subject
        self.fasta = os.path.abspath(fasta)
        self.out = os.path.abspath(out)

        os.makedirs(self.out, exist_ok=True)
        self.script = os.path.join(self.out, "mlst.sh")
        self.log = os.path.join(self.out, "mlst.log")

    def run(self):
        """
        生成并执行
        """
        blastdb = os.path.join(self.database["mlst"], f"pubmlst/{self.subject}/{self.subject}.fa")
        datadir = os.path.join(self.database["mlst"], "pubmlst")
        self.command = f"""# MLST
set -e
source {self.software["activate"]} {self.env["typetrace"]}
mlst --nopath --legacy --quiet --scheme {self.subject} -blastdb {blastdb} -datadir {datadir} \"{self.fasta}\" > {self.out}/mlst.tmp
{self.software["python"]} {BINPATH}/deal_mlst_result.py --input {self.out}/mlst.tmp --scheme {datadir}/{self.subject}/{self.subject}.txt --out {self.out}/mlst.tsv
"""
        with open(self.script, 'w') as OUT:
            print(self.command, file=OUT)
        run(f"{self.software['shell']} {self.script} 2>&1 > {self.log}",
            shell=True)

    def to_upload(self, d_upload):
        """
        link the result file to dir upload
        """
        cmd = f"cp -rf {self.out}/mlst.tsv {d_upload}/mlst.tsv"
        run(cmd, shell=True)

    @property
    def result(self):
        """
        获取mlst分型结果
        """
        f_result = Path(f"{self.out}/mlst.tsv").absolute()
        if f_result.exists():
            with open(f_result, 'r') as IN:
                next(IN)
                for line in IN:
                    arr = line.strip().split("\t")
                    return arr[0]
        else:
            return None
