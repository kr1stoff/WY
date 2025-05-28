#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/9 9:26
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/9 9:26
import logging
import os
import subprocess
from pathlib import Path

import yaml

# 配置文件
d_lib = Path(__file__).absolute().parent
d_base = d_lib.parent
d_config = d_base.joinpath("config")
f_software = d_config.joinpath("software.yml")
software = yaml.safe_load(open(f_software, 'r').read())

# 日志格式
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class PPrimer():
    """
    预扩增引物设计流程
    """

    def __init__(self, reference: Path, genomes: Path, out: Path, run):
        """
        Init the object

        :param reference: The reference sequence
        :param genomes: The other type of the target species
        :param out: The output dir for the result
        :param run: Run the command direct or just generate
        """
        self.d_out = self.set_out(out)
        self.reference = reference.absolute()
        self.genomes = genomes.absolute()
        self.prepare = self.d_out.joinpath("Prepare")
        self.prepare.mkdir(exist_ok=True)
        self.run_direct = run
        self.bin = d_base.joinpath("bin")

        # 常用软件
        self.python = software["python"]
        self.seqkit = software["seqkit"]
        self.bedtools = software["bedtools"]
        self.primer3 = software["primer3"]
        self.mfeprimer = software["mfeprimer"]

        # 软件目录(MicroGMT中需要用到)
        self.microgmt = software["microgmt"]
        self.samtools = software["samtools"]
        self.minimap2 = software["minimap2"]
        self.bcftools = software["bcftools"]

        # 命令
        self.command = []

    def snp(self):
        """
        确认各个突变位点
        """
        self.d_snp = self.d_out.joinpath("SNP")
        self.d_snp.mkdir(exist_ok=True)
        # 环境变量
        cmd = f"""OLDPATH=$PATH
export PATH={self.samtools}:{self.minimap2}:{self.bcftools}:$PATH"""
        self.command.append(cmd)
        # VCF
        d_vcf = self.d_snp.joinpath("VCF")
        d_vcf.mkdir(exist_ok=True)
        # TODO: 添加处理序列名称的步骤
        cmd = f"{self.python} {self.microgmt}/sequence_to_vcf.py -r {self.reference} -i assembly -fs {self.genomes} -o {d_vcf}\nexport PATH=$OLDPATH"
        self.command.append(cmd)

        # Merge
        cmd = f"{self.python} {self.bin}/merge_vcf.py -d {d_vcf} -o {self.d_snp}/all.variants.bed"
        self.command.append(cmd)

    def format(self):
        """
        Upcase the fasta file
        """
        f_out = self.prepare.joinpath("sequence.fasta")
        cmd = f"{self.seqkit} seq {self.reference} -u > {f_out}"
        self.command.append(cmd)
        self.reference = f_out

    def mask(self):
        """
        Mask the useless or snp region
        """
        f_mask = self.d_snp.joinpath("all.variants.bed")
        f_out = self.prepare.joinpath("sequence.final.fa")
        cmd = f"{self.bedtools} maskfasta -soft -fi {self.reference} -bed {f_mask} -fo {f_out}"
        self.reference = f_out
        self.command.append(cmd)

    def internal_design(self, f_config=None):
        """
        内引物及探针的设计
        """
        # TODO: 参数优化
        self.d_internal = self.d_out.joinpath("Internal")
        os.makedirs(self.d_internal, exist_ok=True)

        # 引物预测
        if f_config:
            cmd = f"{self.python} {self.bin}/internal_primer_design.py -f {self.reference} -o {self.d_internal} --config {f_config}"
        else:
            cmd = f"{self.python} {self.bin}/internal_primer_design.py -f {self.reference} -o {self.d_internal}"
        self.command.append(cmd)

        # Primer3 结果整理
        f_out = self.d_internal.joinpath("primers.txt")
        cmd = f"{self.python} {self.bin}/internal_primer_parser.py -f {f_out} --ref {self.reference} -o {self.d_internal}"
        self.command.append(cmd)

    def external_design(self, f_config=None):
        """
        外引物设计

        # TODO: 暂时使用现有方法，后续使用Primer3进行设计
        """
        self.d_external = self.d_out.joinpath("External")
        self.d_external.mkdir(exist_ok=True)

        # 整理内引物的bed文件格式，从而符合primerplex的输入
        cmd = f"{self.python} {self.bin}/internal_to_primerplex.py -i {self.d_internal}/primer.bed -o {self.d_external}/internal.bed"
        self.command.append(cmd)
        # PrimerPlex
        cmd = f"{software['python38']} {software['primerplex']} -ref {self.reference} -region {self.d_external}/internal.bed -th 100 -primernum1 10 -blast -minampllen 175 -maxampllen 250 -optampllen 200 -minprimerlen 18 -maxprimerlen 28 -optprimerlen 23 -minprimermelt 63 -maxprimermelt 67 -optprimermelt 64 -minprimergc 30 -maxprimergc 70 -optprimergc 50 -minprimerendgc 0 -maxprimerendgc 2 -maxprimerpolyn 5 -maxoverlap 300 -returnvariantsnum 5 -skip"
        self.command.append(cmd)
        # 结果整理
        cmd = f"{self.python} {self.bin}/primerplex_parse.py --primerplex {self.d_external}/internal_primers_combination_1_info.xls --ref {self.reference} --prefix {self.d_external}/primer"
        self.command.append(cmd)

    def merge(self):
        """
        合并内外引物设计结果
        """
        self.d_merge = self.d_out.joinpath("Merge")
        self.d_merge.mkdir(exist_ok=True)

        cmd = f"{self.python} {self.bin}/internal_external_merge.py --internal {self.d_internal}/primer.xls --external {self.d_external}/primer.xls --out {self.d_merge}/primer.xls"
        self.command.append(cmd)

    def evaluate(self):
        """
        使用 MFEPrimer 进行引物质量评估
        """
        self.d_mfeprimer = os.path.join(self.d_out, "MEFPrimer")
        os.makedirs(self.d_mfeprimer, exist_ok=True)

        # 为参考建立索引
        cmd = f"{self.mfeprimer} index -i {self.reference} -f"
        self.command.append(cmd)
        # 评估
        # TODO: 参数优化
        f_out = os.path.join(self.d_mfeprimer, "MEFPrimer.result.txt")
        cmd = f"{self.mfeprimer} -d {self.reference} -S 300 -t 55 -i {self.d_primer3}/primer.fasta -o {f_out}"
        self.command.append(cmd)

    def set_out(self, out: Path):
        """
        Set the out put dir

        :param out: The out put dir
        """
        d_out = out.absolute()
        d_out.mkdir(exist_ok=True)
        return out

    def run(self):
        """
        Run the generated script
        """
        logging.info("Start to run")
        command = f"bash {self.f_script} >{self.f_log} 2>&1"
        subprocess.run(command, shell=True)

    def finish(self):
        """
        Finish the pipe
        """
        # 生成脚本
        self.f_script = os.path.join(self.d_out, "run.sh")
        self.f_log = os.path.join(self.d_out, "run.log")
        logging.info(f"Write the command to script {self.f_script}")
        with open(self.f_script, 'w') as OUT:
            print("set -e", file=OUT)
            for command in self.command:
                print(command, file=OUT)

        # 执行脚本
        if self.run_direct:
            self.run()
