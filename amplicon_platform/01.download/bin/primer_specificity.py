#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/17 10:22
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/17 10:22
import logging
import subprocess
import sys
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
class Specificity(object):
    """
    特异性评估
    """

    def __init__(self, f_fasta: Path, f_info: Path, d_out: Path, taxon: str, type_: str,
                 blastn="/usr/bin/blastn", taxonkit="/home/renchaobo/.local/bin/taxonkit",
                 database="/sdbb/bioinfor/yehui/nt/nt.fa", taxondb="/sdbb/bioinfor/yehui/nt", perl="/usr/bin/perl"):
        """
        Init the object
        """
        logger.info(f"Start to generate the command for specificity analysis")
        self.bin = Path(__file__).parent.absolute()
        self.cmd = ["set +e"]
        self.f_fasta = f_fasta
        self.f_info = f_info
        self.d_out = d_out
        self.taxon = taxon
        self.type = type_

        # 软件/数据库
        self.perl = perl
        self.blastn = blastn
        self.db = database
        self.taxondb = taxondb
        self.taxonkit = taxonkit

    def blast(self):
        """
        比对
        """
        self.blastout = self.d_out.joinpath("blast.tsv")
        cmd = f"""# 比对
export BLASTDB={self.taxondb}
{self.blastn} -task blastn-short -query {self.f_fasta} -db {self.db} -num_threads 80 -word_size 9 -evalue 800000 -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs staxid sscinames sstrand' -out {self.blastout}.tmp
{self.perl} {self.bin}/filter_blast.pl {self.blastout}.tmp | sort -k1 - | sed '1i\query\\tseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tqcovs\\ttaxid\\tsscinames\\tsstrand' > {self.blastout} 
"""
        self.cmd.append(cmd)

    def get_taxon_list(self):
        """
        获取包含物种的TaxonID
        """
        self.f_taxonlist = self.d_out.joinpath("taxon.list")
        cmd = f"""# 获取TaxonID列表
{self.taxonkit} list -i {self.taxon} | awk '{{print $1}}' | grep -v '^$' > {self.f_taxonlist}
"""
        self.cmd.append(cmd)

    def single_parse(self):
        """
        单端处理
        """
        self.f_single_stat = self.d_out.joinpath("single.stat")
        cmd = f"""# 单端处理
{self.perl} {self.bin}/judge_primer_blast.pl {self.f_taxonlist} {self.blastout} > {self.f_single_stat}.same
{self.perl} {self.bin}/same_stat.pl {self.f_single_stat}.same > {self.f_single_stat}
"""
        self.cmd.append(cmd)

    def pair_parse(self):
        """
        引物处理
        """
        d_split = self.d_out.joinpath("split")
        self.f_pair_stat = self.d_out.joinpath("pair.stat")
        cmd = f"""# 成对处理
mkdir -p {d_split}
for i in `cut -f 1 {self.f_info}`
do
  {self.perl} {self.bin}/blast_primer_pair.pl {self.f_info} {self.f_single_stat}.same $i > {d_split}/$i.pair.same
  {self.perl} {self.bin}/judge_primer_blast.pl {self.f_taxonlist} {d_split}/$i.pair.same 3 > {d_split}/$i.pair.middle
  {self.perl} {self.bin}/same_stat.pl {d_split}/$i.pair.middle > {d_split}/$i.pair.stat
  {self.perl} {self.bin}/percent_stat.pl $i {self.f_info} {self.f_single_stat} {d_split}/$i.pair.stat > {d_split}/$i.all.stat
done
cat {d_split}/*.all.stat | awk 'NR==1||NR%2==0' > {self.f_pair_stat}
"""
        self.cmd.append(cmd)

    def write_command(self, f_cmd):
        """
        Write the command to the shell script

        :param f_cmd: The output file for the command
        """
        logger.info(f"Write the command to file: {f_cmd}")
        with open(f_cmd, 'w') as OUT:
            print(*self.cmd, sep="\n", file=OUT)

    def final(self, run):
        """

        """
        self.blast()
        self.get_taxon_list()
        self.single_parse()
        if self.type == "primer":
            self.pair_parse()
        self.write_command(self.d_out.joinpath("run.sh"))
        if run:
            command = f"bash {self.d_out.joinpath('run.sh')}"
            subprocess.run(command, shell=True)
            try:
                (pipe_out, pipe_error) = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()
            except:
                sys.exit(logger.error(pipe_error))


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-f', '--fasta',
              required=True,
              type=click.Path(),
              help="The primer/probe fasta file")
@click.option('--info',
              required=True,
              type=click.Path(),
              help="The primer/probe info file")
@click.option('-t', '--type_',
              default="primer",
              show_default=True,
              type=click.Choice(["primer", "probe"]),
              help="The input type")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
@click.option('-d', '--database',
              default="/sdbb/bioinfor/yehui/nt/nt.fa",
              show_default=True,
              type=click.Path(),
              help="The database to use")
@click.option('--taxon',
              required=True,
              help="The target taxon, sep by ','")
@click.option('--blastn',
              default="/usr/bin/blastn",
              show_default=True,
              type=click.Path(),
              help="The blastn path")
@click.option('--taxonkit',
              default="/home/renchaobo/.local/bin/taxonkit",
              show_default=True,
              type=click.Path(),
              help="The taxonkit path")
@click.option('--taxondb',
              default="/sdbb/bioinfor/yehui/nt",
              show_default=True,
              type=click.Path(),
              help="The Path contain taxon info")
@click.option("--run/--no-run",
              default=True,
              show_default=True,
              help="Whether run the pipeline or just generate the script")
def main(fasta, info, type_, out, database, taxon, blastn, taxonkit, taxondb, run):
    """
    特异性评估
    """
    f_fasta = Path(fasta).absolute()
    f_info = Path(info).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)

    pipe = Specificity(f_fasta, f_info, d_out, taxon, type_, blastn=blastn, taxonkit=taxonkit, database=database,
                       taxondb=taxondb)
    pipe.final(run)


if __name__ == "__main__":
    main()
