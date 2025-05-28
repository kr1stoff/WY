#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/11/21 15:53
# @Last Modified by:   Ming
# @Last Modified time: 2023/11/21 15:53
import logging
import subprocess
import sys
from pathlib import Path

import click
import yaml

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Some Setting
PATH = Path(__file__).absolute().parent.parent
CONFIGPATH = PATH.joinpath("config")
BINPATH = PATH.joinpath("bin")
f_software = CONFIGPATH.joinpath("software.yml")
f_database = CONFIGPATH.joinpath("database.yml")
f_env = CONFIGPATH.joinpath("env.yml")


class cgSNP(object):
    """
    微远器械cgSNP分型流程
    
    :TODO add contig support
    """

    def __init__(self, sample, ref, mask, out):
        """
        Init the object

        :param sample: The sample sheet file contain the name and path
        :param ref: The reference fasta file
        :param mask: The mask region bed file
        :param out: The output project
        """
        self.d_out = Path(out).absolute()
        self.d_out.mkdir(exist_ok=True, parents=True)
        self.f_sample = Path(sample).absolute()
        self.f_ref = Path(ref).absolute()
        self.f_mask = Path(mask).absolute()

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

        # 常用设置
        self.parallel = min(len(self.sample), 4)
        self.cpu = 12

        # 软件
        self.software = yaml.safe_load(open(f_software))
        # 环境
        self.env = yaml.safe_load(open(f_env))
        # 数据库
        self.database = yaml.safe_load(open(f_database))

    @property
    def sample(self):
        """
        获取样本名称
        """
        res = []
        with open(self.f_sample, 'r') as IN:
            for line in IN:
                arr = line.strip().split("\t")
                res.append(arr[0])
        return res

    def seqType(self, arr):
        """
        获取序列类型

        :param arr: The line arr
        """
        fasta_suffix = tuple([".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"])
        if len(arr) == 2:
            if arr[1].endswith(fasta_suffix):
                return "fasta"
            else:
                return "se"
        elif len(arr) == 3:
            return "pe"
        else:
            sys.exit(logger.error("Check Your SampleSheet"))

    def prepare(self):
        """
        Prepare the sample info and other info needed follow
        """
        self.d_prepare = self.d_out.joinpath("prepare")
        self.d_prepare.mkdir(exist_ok=True, parents=True)

    def fastp(self):
        """
        Fastp progress

        判断样本序列类型，如果是fasta则直接软链接，如果是fastq，则使用fastp质控
        """
        self.d_fastp = self.d_out.joinpath("fastp")
        self.d_fastp.mkdir(exist_ok=True, parents=True)
        self.f_fastp = self.d_shell.joinpath("fastp.sh")
        self.log_fastp = self.d_log.joinpath("fastp.log")
        with open(self.f_sample, 'r') as IN, open(self.f_fastp, 'w') as OUT:
            for line in IN:
                arr = line.strip().split("\t")
                f_json = f"{self.d_fastp}/{arr[0]}.json"
                f_html = f"{self.d_fastp}/{arr[0]}.html"
                if self.seqType(arr) == "se":
                    cmd = f"{self.software['fastp']} -i {arr[1]} -o {self.d_fastp}/{arr[0]}_1.fq.gz -j {f_json} -h {f_html}"
                elif self.seqType(arr) == "pe":
                    cmd = f"{self.software['fastp']} -i {arr[1]} -I {arr[2]} -o {self.d_fastp}/{arr[0]}_1.fq.gz -O {self.d_fastp}/{arr[0]}_2.fq.gz -j {f_json} -h {f_html}"
                elif self.seqType(arr) == "fasta":
                    cmd = f"ln -sf {arr[1]} {self.d_fastp}/{arr[0]}.fasta"
                else:
                    sys.exit(logger.error("Please Check your SampleSheet"))
                print(cmd, file=OUT)
        self.cmd.append("echo -e Start fastp at `date +%y-%m-%d.%H:%M:%S`")
        self.cmd.append(f"{self.software['parallel']} -j {self.parallel} < {self.f_fastp} >{self.log_fastp} 2>&1")
        self.cmd.append("echo -e End fastp at `date +%y-%m-%d.%H:%M:%S`\n")

    def snippy(self):
        """
        Use Snippy to call SNP
        """
        self.d_snippy = self.d_out.joinpath("snippy")
        self.d_snippy.mkdir(exist_ok=True, parents=True)
        self.f_snippy = self.d_shell.joinpath("snippy.sh")
        self.log_snippy = self.d_log.joinpath("snippy.log")
        with open(self.f_sample, 'r') as IN, open(self.f_snippy, 'w') as OUT:
            for line in IN:
                arr = line.strip().split("\t")
                sample = arr[0]
                d_out = self.d_snippy.joinpath("sample")
                if self.seqType(arr) == "pe":
                    r1 = self.d_fastp.joinpath(f"{sample}_1.fq.gz")
                    r2 = self.d_fastp.joinpath(f"{sample}_2.fq.gz")
                    cmd = f"source {self.software['activate']} {self.env['snippy']};snippy --force --cpus {self.cpu} --ref {self.f_ref} --R1 {r1} --R2 {r2} --outdir {self.d_snippy}/{sample} --prefix {sample}"
                elif self.seqType(arr) == "se":
                    r1 = self.d_fastp.joinpath(f"{sample}_1.fq.gz")
                    cmd = f"source {self.software['activate']} {self.env['snippy']};snippy --force --cpus {self.cpu} --ref {self.f_ref} --se {r1} --outdir {self.d_snippy}/{sample} --prefix {sample}"
                elif self.seqType(arr) == "fasta":
                    ctg = self.d_fastp.joinpath(f"{sample}.fasta")
                    cmd = f"source {self.software['activate']} {self.env['snippy']};snippy --force --cpus {self.cpu} --ref {self.f_ref} --ctgs {ctg} --outdir {self.d_snippy}/{sample} --prefix {sample}"
                else:
                    sys.exit(logger.error("Check Your SampleSheet"))

                print(cmd, file=OUT)
        self.cmd.append("echo -e Start snippy at `date +%y-%m-%d.%H:%M:%S`")
        self.cmd.append(f"{self.software['parallel']} -j {self.parallel} < {self.f_snippy} >{self.log_snippy} 2>&1")
        self.cmd.append("echo -e End snippy at `date +%y-%m-%d.%H:%M:%S`\n")

    def post_snippy(self):
        """
        Merge the snippy result
        """
        self.f_postsnippy = self.d_shell.joinpath("post_snippy.sh")
        self.log_postsnippy = self.d_log.joinpath("post_snippy.log")
        with open(self.f_postsnippy, 'w') as OUT:
            cmd = f"""cd {self.d_snippy}
source {self.software['activate']} {self.env['snippy']}
snippy-core --ref {self.f_ref} {' '.join(self.sample)} --mask {self.f_mask}
snippy-clean_full_aln core.full.aln > gubbins.aln
"""
            print(cmd, file=OUT)
        self.cmd.append("echo -e Start post-snippy at `date +%y-%m-%d.%H:%M:%S`")
        self.cmd.append(f"bash {self.f_postsnippy} >{self.log_postsnippy} 2>&1")
        self.cmd.append("echo -e End post-snippy at `date +%y-%m-%d.%H:%M:%S`\n")

    def gubbins(self):
        """
        Use gubbins to regerate the core region of snp
        """
        self.d_gubbins = self.d_out.joinpath("gubbins")
        self.d_gubbins.mkdir(exist_ok=True, parents=True)
        self.f_gubbins = self.d_shell.joinpath("gubbins.sh")
        self.log_gubbins = self.d_log.joinpath("gubbins.log")
        with open(self.f_gubbins, 'w') as OUT:
            cmd = f"""cd {self.d_gubbins}
source {self.software['activate']} {self.env['gubbins']}
run_gubbins.py --filter-percentage 50 --threads {self.cpu} {self.d_snippy}/gubbins.aln
"""
            print(cmd, file=OUT)
        self.cmd.append("echo -e Start gubbins at `date +%y-%m-%d.%H:%M:%S`")
        self.cmd.append(f"bash {self.f_gubbins} >{self.log_gubbins} 2>&1")
        self.cmd.append("echo -e End gubbins at `date +%y-%m-%d.%H:%M:%S`\n")

    def plot(self):
        """
        Plot the tree
        """
        self.d_plot = self.d_out.joinpath("plot")
        self.d_plot.mkdir(exist_ok=True, parents=True)
        self.f_plot = self.d_shell.joinpath("plot.sh")
        self.log_plot = self.d_log.joinpath("plot.log")
        with open(self.f_plot, 'w') as OUT:
            cmd = f"""cd {self.d_plot}
{self.software['rscript']} {BINPATH}/tree.R --branchlength --tree {self.d_gubbins}/gubbins.final_tree.tre --prefix cgSNP
"""
            print(cmd, file=OUT)
        self.cmd.append("echo -e Start plot at `date +%y-%m-%d.%H:%M:%S`")
        self.cmd.append(f"bash {self.f_plot} >{self.log_plot} 2>&1")
        self.cmd.append("echo -e End plot at `date +%y-%m-%d.%H:%M:%S`\n")

    def finish(self, run: bool = False):
        """

        """
        logger.info(f"Write the command to {self.f_shell}")
        with open(self.f_shell, 'w') as OUT:
            for cmd in self.cmd:
                print(cmd, file=OUT)

        if run:
            logger.info(f"Start to run the pipe")
            handle_log = open(self.f_log, 'w')
            result = subprocess.run(["bash", str(self.f_shell)], stdout=handle_log, stderr=handle_log)
            handle_log.close()
            if result.returncode == 0:
                return 0
            else:
                return -1

    def run(self):
        """
        Run the pipe
        """
        logger.info(f"Start to run the pipe")
        self.prepare()
        self.fastp()
        self.snippy()
        self.post_snippy()
        self.gubbins()
        self.plot()
        self.finish(run=True)

    def to_upload(self, d_upload):
        """
        copy the result file to dir upload
        """
        cmd = f"""set -e
cp -rf {self.d_out}/gubbins/gubbins.final_tree.tre {d_upload}/cgSNP.tre
cp -rf {self.d_out}/plot/cgSNP.png {d_upload}/cgSNP.png
"""
        subprocess.run(cmd, shell=True, check=True)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-s", "--sample",
              required=True,
              type=click.Path(),
              help="The sample sheet file")
@click.option("-r", "--reference",
              required=True,
              type=click.Path(),
              help="The reference file")
@click.option("-m", "--mask",
              required=True,
              type=click.Path(),
              help="The mask region bed file")
@click.option("-o", "--out",
              type=click.Path(),
              help="The out put dir")
@click.option("--run/--no-run",
              default=True,
              show_default=True,
              help="Whether run the pipe direct")
def cli(sample, reference, mask, out, run):
    """
    cgSNP流程
    """
    f_sample = Path(sample).absolute()
    f_reference = Path(reference).absolute()
    f_mask = Path(mask).absolute()
    d_out = Path(out).absolute()

    project = cgSNP(f_sample, f_reference, f_mask, d_out)
    project.prepare()
    project.fastp()
    project.snippy()
    project.post_snippy()
    project.gubbins()
    project.plot()
    project.finish(run=run)


if __name__ == "__main__":
    cli()
