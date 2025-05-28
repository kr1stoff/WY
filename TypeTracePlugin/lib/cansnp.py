#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2024/1/9 15:27
# @Last Modified by:   Ming
# @Last Modified time: 2024/1/9 15:27
import logging
import subprocess
import sys
from pathlib import Path

import click
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


class CanSNP(object):
    """
    CanSNP pipeline
    """

    def __init__(self, fasta, reference, database, out):
        """
        Init the object

        :param fasta: 组装得到的基因组序列
        :param reference: 参考所在目录
        """
        self.f_fasta = Path(fasta).absolute()
        self.d_out = Path(out).absolute()
        self.d_out.mkdir(exist_ok=True, parents=True)
        self.d_reference = Path(reference).absolute()
        self.f_database = Path(database).absolute()

        # 常用目录
        self.d_log = self.d_out.joinpath("log")
        self.d_log.mkdir(exist_ok=True, parents=True)
        self.d_shell = self.d_out.joinpath("shell")
        self.d_shell.mkdir(exist_ok=True, parents=True)
        # 总脚本
        self.cmd = []
        self.cmd.append("set -e")
        self.f_shell = self.d_shell.joinpath("all.sh")
        self.f_log = self.d_log.joinpath("all.log")

        # 软件
        self.software = yaml.safe_load(open(f_software))
        # 数据库
        self.env = yaml.safe_load(open(f_env))
        # 环境
        self.database = yaml.safe_load(open(f_database))

    def prepare(self):
        """
        Link the database and sequence
        """
        self.d_prepare = self.d_out.joinpath("prepare")
        self.d_prepare.mkdir(exist_ok=True, parents=True)
        self.f_prepare = self.d_shell.joinpath("prepare.sh")
        self.log_prepare = self.d_log.joinpath("prepare.log")
        cmd = f"""set -e
ln -sf {self.f_fasta} {self.d_prepare}
"""
        with open(self.f_prepare, 'w') as OUT:
            print(cmd, file=OUT)
        self.cmd.append("echo -e Start prepare at `date +%y-%m-%d.%H:%M:%S`")
        self.cmd.append(f"bash {self.f_prepare} >{self.log_prepare} 2>&1")
        self.cmd.append("echo -e End perpare at `date +%y-%m-%d.%H:%M:%S`\n")

    def cansnp(self):
        """
        Run the CanSNP analysis with software cansnper2
        """
        self.d_cansnp = self.d_out.joinpath("cansnp")
        self.d_cansnp.mkdir(exist_ok=True, parents=True)
        self.f_cansnp = self.d_shell.joinpath("cansnp.sh")
        self.log_cansnp = self.d_log.joinpath("cansnp.log")
        cmd = f"""set -e
cd {self.d_cansnp}
source {self.software['activate']} {self.env['cansnper2']}
CanSNPer2 --database {self.f_database} \
--refdir {self.d_reference} \
--workdir {self.d_cansnp} \
--tmpdir {self.d_cansnp}/tmp \
{self.d_prepare}/{self.f_fasta.name}
awk -F' ' -vOFS=\"\\t\" 'BEGIN{{print "Name","CanSNP"}}$1=="Final"{{print "Query",$3}}' {self.d_cansnp}/results/*_snps.txt > {self.d_cansnp}/cansnp.tsv
"""
        with open(self.f_cansnp, 'w') as OUT:
            print(cmd, file=OUT)
        self.cmd.append("echo -e Start CanSNP at `date +%y-%m-%d.%H:%M:%S`")
        self.cmd.append(f"bash {self.f_cansnp} >{self.log_cansnp} 2>&1")
        self.cmd.append("echo -e End CanSNP at `date +%y-%m-%d.%H:%M:%S`\n")

    def finish(self, run: bool = False):
        """
        Write the command to the script and exec
        """
        logger.info(f"Write the command to {self.f_shell}")
        with open(self.f_shell, 'w') as OUT:
            for cmd in self.cmd:
                print(cmd, file=OUT)

        if run:
            logger.info(f"Start to run the pipe")
            try:
                subprocess.run(f"bash {self.f_shell}", shell=True, check=True)
            except Exception as e:
                logger.error(e)
                sys.exit(1)

    def run(self):
        """
        执行
        """
        self.prepare()
        self.cansnp()
        self.finish(run=True)

    def to_upload(self, d_upload):
        """
        copy the result file to dir upload
        """
        cmd = f"""set -e
cp -rf {self.d_out}/cansnp/results/*snps.txt {d_upload}/cansnp_detail.txt
cp -rf {self.d_out}/cansnp/cansnp.tsv {d_upload}/cansnp.tsv
"""
        subprocess.run(cmd, shell=True, check=True)

    @property
    def result(self):
        """
        获取canSNP分型结果
        """
        f_result = Path(f"{self.d_out}/cansnp/cansnp.tsv").absolute()
        if f_result.exists():
            with open(f_result, 'r') as IN:
                next(IN)
                for line in IN:
                    arr = line.strip().split("\t")
                    return arr[1]
        else:
            return None


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-f', '--fasta',
              required=True,
              type=click.Path(),
              help="The input fasta genome file")
@click.option('-db', '--database',
              required=True,
              type=click.Path(),
              help="The database to use")
@click.option('-r', '--ref',
              required=True,
              type=click.Path(),
              help="The reference dir")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
@click.option("--run/--no-run",
              default=True,
              show_default=True,
              help="Whether run the pipe direct")
def main(fasta, database, ref, out, run):
    """
    CanSNP 流程
    """
    f_fasta = Path(fasta).absolute()
    f_db = Path(database).absolute()
    d_ref = Path(ref).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(parents=True, exist_ok=True)

    job = CanSNP(f_fasta, d_ref, f_db, d_out)
    job.prepare()
    job.cansnp()
    job.finish(run)


if __name__ == "__main__":
    main()
