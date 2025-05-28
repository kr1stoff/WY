#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/12/8 10:58
# @Last Modified by:   Ming
# @Last Modified time: 2022/12/8 10:58
import logging
import subprocess
import sys
from pathlib import Path

import click
import yaml
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

#### Some Global Variable
f_software = Path(__file__).absolute().parent.parent.joinpath("config/software.yml")
software = yaml.safe_load(open(f_software, 'r'))
bin = Path(__file__).absolute().parent


#### Functions
class Pipe(object):
    """
    病毒保守结构域预测
    """

    def __init__(self, name: str, ref: Path, f_genome: Path, d_out: Path, cpu: int, methods: list, merge: int,
                 merge_max: int):
        """
        Init the object

        :param name: The latin name
        :param ref: The referece genome
        :param f_genome: The file contain the genomes
        :param d_out: The output dir
        :param cpu: The cpu number to use
        """
        self.name = name
        self.f_ref = ref.absolute()
        self.f_genome = f_genome.absolute()
        self.d_out = d_out.absolute()
        self.d_out.mkdir(exist_ok=True, parents=True)
        self.cpu = cpu
        self.methods = methods
        self._merge = merge
        self._merge_max = merge_max

    def prepare(self):
        """
        文件预处理
        """
        self.d_prepare = self.d_out.joinpath("prepare")
        self.d_prepare.mkdir(exist_ok=True, parents=True)

        logger.info(f"Merge the reference and genomes")
        self.f_ref_format = self.d_prepare.joinpath("ref.fna")
        self.f_merge = self.d_prepare.joinpath("merge.fna")
        self.f_merge_noref = self.d_prepare.joinpath("merge.noref.fna")
        with open(self.f_merge, 'w') as OUT1, open(self.f_merge_noref, 'w') as OUT2, open(self.f_ref_format,
                                                                                          'w') as OUT3:
            for record in SeqIO.parse(self.f_ref, 'fasta'):
                print(f">{record.name}\n{record.seq}", file=OUT1)
                print(f">{record.name}\n{record.seq}", file=OUT3)
                self.ref_name = record.name
            for record in SeqIO.parse(self.f_genome, "fasta"):
                if record.name != self.ref_name:
                    print(f">{record.name}\n{record.seq}", file=OUT1)
                    print(f">{record.name}\n{record.seq}", file=OUT2)

    def phast_pipe(self):
        """
        Use PHAST to predict conserved region

        TODO: 软件参数需要进一步评估
        """
        logger.info(f"PHAST Analysis")
        self.d_phast = self.d_out.joinpath("phast")
        self.d_phast.mkdir(exist_ok=True)

        cmd = f"""set -e
{software['mafft']} --quiet --thread {self.cpu} --6merpair --keeplength --add {self.f_merge_noref} {self.f_ref_format} > {self.d_phast}/align.fa
{software['rapidnj']} -c 24 -a jc -t d -i fa {self.d_phast}/align.fa -b 100 -x {self.d_phast}/align.nwk
sed -i "s/'//g" {self.d_phast}/align.nwk
{software['phyloFit']} --quiet --tree {self.d_phast}/align.nwk --subst-mod HKY85 --out-root {self.d_phast}/modelfile {self.d_phast}/align.fa
{software['phastCons']} --quiet -N {self.ref_name} --most-conserved {self.d_phast}/mostcons.bed {self.d_phast}/align.fa {self.d_phast}/modelfile.mod > {self.d_phast}/scores.wig
"""
        logger.debug(cmd)

        try:
            (phast_out, phast_error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
        except:
            sys.exit(logger.error(phast_error))

        # Merge neighbor
        cmd = f"{software['python']} {bin}/deal_phast_bed.py --phast {self.d_phast} --chromosomes {self.ref_name} --name '{self.name}' --merge {self._merge} --merge_max {self._merge_max} --out {self.d_phast}/phast.bed"
        logger.debug(cmd)
        try:
            (merge_out, merge_error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
        except:
            sys.exit(logger.error(merge_error))

        subprocess.run(cmd, shell=True)

    def pairedkmers_pipe(self):
        """
        Use pared-kmers to predict conserved region
        """
        logger.info(f"paired-kmers Analysis")
        self.d_pairedkmers = self.d_out.joinpath("pairedkmers")
        self.d_pairedkmers.mkdir(exist_ok=True)

        cmd = f"{software['paired-kmers']} -i {self.f_merge} -k 17 -s 90 -S 300 -c 24 -o {self.d_phast}/kmer.bed"
        logger.debug(cmd)
        try:
            (pairedkmers_out, pairedkmers_error) = subprocess.Popen(cmd, shell=True,
                                                                    stdout=subprocess.PIPE).communicate()
        except:
            sys.exit(logger.error(pairedkmers_error))
        subprocess.run(cmd, shell=True)

        cmd = f"{software['python']} {bin}/deal_pairedkmers_bed.py --pairedkmers {self.d_pairedkmers} --chromosomes {self.ref_name} --name '{self.name}' --out {self.d_pairedkmers}/pkmers.bed"
        logger.debug(cmd)
        try:
            (merge_out, merge_error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()
        except:
            sys.exit(logger.error(merge_error))
        subprocess.run(cmd, shell=True)

    def run(self):
        """
        流程执行
        """
        self.prepare()
        if "phast" in self.methods:
            self.phast_pipe()
        if "pairedkmers" in self.methods:
            self.pairedkmers_pipe()


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.1")
@click.option('--name',
              required=True,
              help="The Latin name")
@click.option('--ref',
              required=True,
              type=click.Path(),
              help="The viral genome file used as reference")
@click.option('--genome',
              required=True,
              type=click.Path(),
              help="The all ref genome fasta file(需保证一个基因组一个文件)")
@click.option("-o", "--out",
              required=True,
              type=click.Path(),
              help="The out put dir")
@click.option("--methods",
              default="phast",
              type=click.Choice(["phast", "pairedkmers", "all"]),
              help="The method used to predict conserved region")
@click.option("--cpu",
              default=64,
              type=int,
              help="The cpu number to use")
@click.option('--merge',
              default=50,
              show_default=True,
              type=int,
              help="当区间距离小于此值时合并区间")
@click.option('--merge_max',
              default=600,
              show_default=True,
              type=int,
              help="所合并区间的最大长度")
def main(name, ref, genome, out, methods, cpu, merge, merge_max):
    """
    病毒保守区域寻找
    """
    f_ref = Path(ref).absolute()
    f_genome = Path(genome).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)
    if methods == "all":
        methods = ["phast", "pairedkmers"]
    else:
        methods = [methods]

    pipe = Pipe(name, f_ref, f_genome, d_out, cpu, methods, merge, merge_max)
    pipe.run()


if __name__ == "__main__":
    main()
