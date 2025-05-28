#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2024/2/22 13:56
# @Last Modified by:   Ming
# @Last Modified time: 2024/2/22 13:56
import logging
import subprocess
from pathlib import Path

import yaml

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Some Setting
PATH = Path(__file__).absolute().parent.parent
CONFIGPATH = PATH.joinpath("config")
BINPATH = PATH.joinpath("bin")
f_software = CONFIGPATH.joinpath("software.yml")
f_database = CONFIGPATH.joinpath("database.yml")
f_env = CONFIGPATH.joinpath("env.yml")


class MLVA(object):
    """
    MLVA 分析流程
    """

    def __init__(self, fasta: Path, primer: Path, out: Path):
        """
        Init the object

        :param fasta: The input genome fasta file
        :param primer: The primer file
        :param out: The output dir
        """
        self.env = yaml.safe_load(open(f_env, 'r').read())
        self.software = yaml.safe_load(open(f_software, 'r').read())
        self.database = yaml.safe_load(open(f_database, 'r').read())

        self.fasta = fasta.absolute()
        self.primer = primer.absolute()
        self.d_out = out.absolute()
        self.d_out.mkdir(parents=True, exist_ok=True)
        self.script = self.d_out.joinpath("mlva.sh")
        self.log = self.d_out.joinpath("mlva.log")

    def run(self):
        """
        生成脚本并执行
        """
        self.command = f"""# MLVA
set -e
mkdir -p {self.d_out}/mlva
ln -sf {self.fasta} {self.d_out}/mlva/Query.fa
source {self.software["activate"]} {self.env["mlva"]}
python {BINPATH}/mlva.py -i {self.d_out}/mlva -o {self.d_out} -p {self.primer}
cat {self.d_out}/MLVA_analysis_mlva.csv | tr ',' '\t' > {self.d_out}/mlva.tsv
"""
        with open(self.script, 'w') as OUT:
            print(self.command, file=OUT)
        subprocess.run(f"{self.software['shell']} {self.script} > {self.log} 2>&1",
                       shell=True, check=True)

    def to_upload(self, d_upload):
        """
        copy the result file to dir upload
        """
        cmd = f"""set -e
cp -rf {self.d_out}/mlva.tsv {d_upload}/
"""
        subprocess.run(cmd, shell=True, check=True)
