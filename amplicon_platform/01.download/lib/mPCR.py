#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/9 9:26
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/9 9:26
import logging
import os
import subprocess

import yaml

# 配置文件
BIN = os.path.dirname(os.path.abspath(__file__))
d_config = os.path.abspath(f"{BIN}/../config")
f_software = os.path.join(d_config, "software.yml")
software = yaml.safe_load(open(f_software, 'r').read())

# 日志格式
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class Primer():
    """
    微远器械多重PCR引物设计流程
    """

    def __init__(self, reference, genomes, out, run):
        """
        Init the object

        :param reference: The reference sequence
        :param genomes: The other type of the target species
        :param out: The out put dir for the result
        :param run: Run the command direct or just generate
        """
        self.d_out = self.set_out(out)
        self.reference = os.path.abspath(reference)
        self.genomes = os.path.abspath(genomes)
        self.prepare = os.path.join(self.d_out, "Prepare")
        os.makedirs(self.prepare, exist_ok=True)
        self.run_direct = run
        self.bin = os.path.join(BIN, "../bin")

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
        d_out = os.path.join(self.d_out, "SNP")
        self.d_snp = d_out
        os.makedirs(d_out, exist_ok=True)
        # 环境变量
        cmd = f"""OLDPATH=$PATH
export PATH={self.samtools}:{self.minimap2}:{self.bcftools}:$PATH"""
        self.command.append(cmd)
        # VCF
        d_vcf = os.path.join(d_out, "VCF")
        os.makedirs(d_vcf, exist_ok=True)
        # TODO: 添加处理序列名称的步骤
        cmd = f"{self.python} {self.microgmt}/sequence_to_vcf.py -r {self.reference} -i assembly -fs {self.genomes} -o {d_vcf}\nexport PATH=$OLDPATH"
        self.command.append(cmd)

        # Merge
        cmd = f"for i in {d_vcf}/*.vcf;do grep -v '#' $i;done > {d_out}/all.variants.vcf"
        self.command.append(cmd)

        # Convert to bed
        cmd = f"awk '{{if(length($4)>length($5))print $1\"\\t\"($2-1)\"\\t\"($2+length($4)-1);else print $1\"\\t\"($2-1)\"\\t\"($2+length($5)-1)}}' {d_out}/all.variants.vcf | {self.bedtools} sort | {self.bedtools} merge > {d_out}/all.variants.bed"
        self.command.append(cmd)

    def format(self):
        """
        Upcase the fasta file
        """
        f_out = os.path.join(self.prepare, "sequence.fasta")
        cmd = f"{self.seqkit} seq {self.reference} -u > {f_out}"
        self.command.append(cmd)
        self.reference = f_out

    def mask(self):
        """
        Mask the useless or snp region
        """
        f_mask = os.path.join(self.d_snp, "all.variants.bed")
        f_out = os.path.join(self.prepare, "sequence.final.fa")
        cmd = f"{self.bedtools} maskfasta -soft -fi {self.reference} -bed {f_mask} -fo {f_out}"
        self.reference = f_out
        self.command.append(cmd)

    def design(self, f_config=None):
        """
        使用 Primer3 进行引物设计

        :param f_config: The project YAML config file
        """
        # TODO: 参数优化
        self.d_primer3 = os.path.join(self.d_out, "Primer3")
        os.makedirs(self.d_primer3, exist_ok=True)

        # 引物预测
        if f_config:
            cmd = f"{self.python} {self.bin}/p3_run.py -f {self.reference} -o {self.d_primer3} --config {f_config}"
        else:
            cmd = f"{self.python} {self.bin}/p3_run.py -f {self.reference} -o {self.d_primer3}"
        self.command.append(cmd)

        # Primer3 结果整理
        f_out = os.path.join(self.d_primer3, "primers.txt")
        cmd = f"{self.python} {self.bin}/p3_parser.py -f {f_out} -o {self.d_primer3}"
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

    def set_out(self, out):
        """
        Set the out put dir

        :param out: The out put dir
        """
        out = os.path.abspath(out)
        os.makedirs(out, exist_ok=True)
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
