#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/8/7 10:07
# @Last Modified by:   Ming
# @Last Modified time: 2023/8/7 10:07
import logging
import subprocess
import sys
from pathlib import Path

import yaml

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# 配置文件
d_lib = Path(__file__).absolute().parent
d_base = d_lib.parent
d_config = d_base.joinpath("config")
f_software = d_config.joinpath("software.yml")
software = yaml.safe_load(open(f_software, 'r').read())


class Sequence(object):
    """
    测序引物设计
    """

    def __init__(self, reference: Path, genomes: Path, out: Path, run):
        """
        Init the object

        :param reference: The reference sequence
        :param genomes: The other type of the target species
        :param out: The output dir for the result
        :param run: Run the command direct or just generate
        """
        self.d_out = Path(out).absolute()
        self.d_out.mkdir(exist_ok=True, parents=True)
        self.reference = Path(reference).absolute()
        self.genomes = Path(genomes).absolute()
        self.bin = d_base.joinpath("bin")
        self.run_direct = run

        # 常用软件
        self.python = software["python"]
        self.perl = software["perl"]
        self.seqkit = software["seqkit"]
        self.bedtools = software["bedtools"]
        self.primer3 = software["primer3"]
        self.mfeprimer = software["mfeprimer"]
        self.blastn = software["blastn"]

        # 软件目录(MicroGMT中需要用到)
        self.microgmt = software["microgmt"]
        self.samtools = software["samtools"]
        self.minimap2 = software["minimap2"]
        self.bcftools = software["bcftools"]

        # 软件目录(MicroGMT中需要用到)
        self.microgmt = software["microgmt"]
        self.samtools = software["samtools"]
        self.minimap2 = software["minimap2"]
        self.bcftools = software["bcftools"]

        # 命令
        self.command = ["set -e"]

    def prepare(self):
        """
        Upcase the fasta file
        """
        self.d_prepare = self.d_out.joinpath("prepare")
        self.d_prepare.mkdir(exist_ok=True, parents=True)

        # 参考处理
        f_out = self.d_prepare.joinpath("sequence.fasta")
        cmd = f"""echo -e [`date +%y-%m-%d.%H:%M:%S`] Start Prepare
{self.seqkit} seq {self.reference} -u > {f_out}"""
        self.command.append(cmd)
        self.reference = f_out
        self.ori = f_out

        # 基因组序列处理
        f_genomes_rename = self.d_prepare.joinpath("genomes.fasta")
        cmd = f"{self.python} {self.bin}/genome_rename.py --fasta {self.genomes} --out {f_genomes_rename}"
        self.command.append(cmd)
        self.command.append("echo -e [`date +%y-%m-%d.%H:%M:%S`] End Prepare\n")
        self.genomes_rename = f_genomes_rename

    def snp(self, freq: float = 0.05):
        """
        确认各个突变位点

        :param freq: The freq of a snp to be regard as a real site
        """
        self.d_snp = self.d_out.joinpath("SNP")
        self.d_snp.mkdir(exist_ok=True)
        # 环境变量
        cmd = f"""echo -e [`date +%y-%m-%d.%H:%M:%S`] Start SNP
OLDPATH=$PATH
export PATH={self.samtools}:{self.minimap2}:{self.bcftools}:$PATH"""
        self.command.append(cmd)
        # VCF
        d_vcf = self.d_snp.joinpath("VCF")
        d_vcf.mkdir(exist_ok=True)
        # TODO: 添加处理序列名称的步骤
        cmd = f"{self.python} {self.microgmt}/sequence_to_vcf.py -t 40 -r {self.reference} -i assembly -fs {self.genomes_rename} -o {d_vcf}\nexport PATH=$OLDPATH"
        self.command.append(cmd)

        # Merge
        cmd = f"{self.python} {self.bin}/merge_vcf.py -f {freq} -d {d_vcf} -o {self.d_snp}/all.variants.bed"
        self.command.append(cmd)

        # Mask
        f_mask = self.d_snp.joinpath("all.variants.bed")
        f_out = self.d_prepare.joinpath("sequence.final.fa")
        cmd = f"{self.bedtools} maskfasta -soft -fi {self.reference} -bed {f_mask} -fo {f_out}"
        self.reference = f_out
        self.command.append(cmd)
        self.command.append("echo -e [`date +%y-%m-%d.%H:%M:%S`] End SNP\n")

    def design(self, region: Path, config: Path, freq: float = 0.05):
        """
        测序引物设计

        :param region: The region bed file of your designed primer amplified region
        :param config: The config file of Sequence
        """
        self.d_design = self.d_out.joinpath("Design")
        self.d_design.mkdir(exist_ok=True, parents=True)
        self.f_region = region.absolute()
        self.f_config = config.absolute()

        # design
        cmd = f"""echo -e [`date +%y-%m-%d.%H:%M:%S`] Start Design
{self.python} {self.bin}/sequence_primer_design.py -r {self.f_region} -f {self.reference} -o {self.d_design} --config {self.f_config}"""
        self.command.append(cmd)

        # degenerate
        cmd = f"{self.python} {self.bin}/sequence_primer_degenerate.py --merge {self.d_design}/primer.xls --snp {self.d_snp}/all.variants.bed -o {self.d_design}/primer.degenerate.xls --freq {freq}"
        self.command.append(cmd)

        # info
        cmd = f"{self.python} {self.bin}/primer_split.py --data {self.d_design}/primer.xls -o {self.d_design}"
        self.command.append(cmd)
        self.command.append("echo -e [`date +%y-%m-%d.%H:%M:%S`] END Design\n")

    def specificity(self, taxonid):
        """
        特异性分析
        """
        if taxonid == "None":
            self.d_specificity = None
        else:
            self.d_specificity = self.d_out.joinpath("Specificity")
            self.d_specificity.mkdir(exist_ok=True, parents=True)
            cmd = f"""echo -e [`date +%y-%m-%d.%H:%M:%S`] Start Specificity
{self.python} {self.bin}/primer_specificity.py -f {self.d_design}/sequence.fasta --info {self.d_design}/primer.info.tsv -t primer -o {self.d_specificity} --taxon {taxonid}
echo -e [`date +%y-%m-%d.%H:%M:%S`] End Specificity
"""
            self.command.append(cmd)

    def inclusive(self, database):
        """
        包容性分析
        """
        self.d_inclusive = self.d_out.joinpath("Inclusive")
        self.d_inclusive.mkdir(exist_ok=True, parents=True)

        f_chromosome_info = Path(database).absolute().parent.joinpath("assembly_chromosome.tsv")
        cmd = f"""echo -e [`date +%y-%m-%d.%H:%M:%S`] Start Inclusive
{self.blastn} -task blastn-short -query {self.d_design}/sequence.fasta -db {database} -num_threads 60 -word_size 9 -evalue 800000 -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sstrand' -out {self.d_inclusive}/blast.tmp
{self.perl} {self.bin}/filter_blast.pl {self.d_inclusive}/blast.tmp > {self.d_inclusive}/blast.tsv
{self.python} {self.bin}/primer_inclusive.py --blast {self.d_inclusive}/blast.tsv --info {self.d_design}/primer.info.tsv -t primer --chromosome_info {f_chromosome_info} --multi -o {self.d_inclusive}/inclusive.tsv
echo -e [`date +%y-%m-%d.%H:%M:%S`] End Inclusive
"""
        self.command.append(cmd)

    def merge(self, taxonname):
        """
        Merge the result
        """
        self.d_merge = self.d_out.joinpath("Merge")
        self.d_merge.mkdir(exist_ok=True, parents=True)

        param = "" if self.d_specificity is None else f"-s {self.d_specificity}/pair.stat"

        cmd = f"""# Merge the result
echo -e [`date +%y-%m-%d.%H:%M:%S`] Start Merge
{self.python} {self.bin}/sequence_primer_merge.py --info {self.d_design}/primer.degenerate.xls -i {self.d_inclusive}/inclusive.tsv {param} --species {taxonname} -o {self.d_merge}/merge.tsv
{self.python} {self.bin}/sequence_primer_to_other_format.py -i {self.d_merge}/merge.tsv -o {self.d_merge}/{taxonname}.sequence.check.tsv -r {self.ori}
echo -e [`date +%y-%m-%d.%H:%M:%S`] End Merge
"""
        self.command.append(cmd)

    def run(self):
        """
        Run the generated script
        """
        logger.info("Start to run")
        command = f"bash {self.f_script} > {self.f_log} 2>{self.f_err}"
        subprocess.run(command, shell=True)
        try:
            (pipe_out, pipe_error) = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()
        except:
            sys.exit(logger.error(pipe_error))

    def finish(self):
        """
        Finish the pipe
        """
        # 生成脚本
        self.f_script = self.d_out.joinpath("run.sh")
        self.f_log = self.d_out.joinpath("run.log")
        self.f_err = self.d_out.joinpath("run.err")
        logger.info(f"Write the command to script {self.f_script}")
        with open(self.f_script, 'w') as OUT:
            print("set -e", file=OUT)
            for command in self.command:
                print(command, file=OUT)

        # 执行脚本
        if self.run_direct:
            self.run()
