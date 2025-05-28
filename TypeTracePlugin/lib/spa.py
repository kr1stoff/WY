#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2024/1/4 9:35
# @Last Modified by:   Ming
# @Last Modified time: 2024/1/4 9:35
import logging
import subprocess
import sys
from pathlib import Path

import click
import yaml
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Some Setting
PATH = Path(__file__).absolute().parent
CONFIGPATH = PATH.parent.joinpath("config")
f_software = CONFIGPATH.joinpath("software.yml")
f_database = CONFIGPATH.joinpath("database.yml")
f_env = CONFIGPATH.joinpath("env.yml")


class Spa(object):
    """
    金黄色葡萄球菌的Spa分型流程
    """

    def __init__(self, fasta, f_primer, db_spa, type_spa, d_out, mismatches=2, minlen=100, maxlen=1000):
        """
        Init the object

        :param fasta: The input fasta genome file
        :param f_primer: The input primer fasta file
        :param db_spa: The
        :param type_spa: The
        :param d_out: The output directory
        :param mismatches: The mismatch count allowed per primer
        :param f_primer: The input primer fasta file
        :param minlen: The minimum required selection length
        :param maxlen: The maximum allowed selection length
        """
        self.fasta = Path(fasta).absolute()
        self.d_out = Path(d_out).absolute()
        self.d_out.mkdir(exist_ok=True, parents=True)
        self.f_primer = Path(f_primer).absolute()
        self.db_spa = Path(db_spa).absolute()
        self.type_spa = Path(type_spa).absolute()
        self.mismatches = mismatches
        self.minlen = minlen
        self.maxlen = maxlen
        self.f_out = self.d_out.joinpath("spa.tsv")

        self.env = yaml.safe_load(open(f_env, 'r').read())
        self.software = yaml.safe_load(open(f_software, 'r').read())
        self.database = yaml.safe_load(open(f_database, 'r').read())

    def prepare(self):
        """
        Prepare the pipe
        """
        self.d_prepare = self.d_out.joinpath("prepare")
        self.d_prepare.mkdir(exist_ok=True, parents=True)
        self.f_fasta = self.d_prepare.joinpath("Query.fasta")
        cmd = f"ln -sf {self.fasta} {self.f_fasta}"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        ## make blast index
        logger.info(f"make index for {self.f_fasta}")
        blastdb = self.d_prepare.joinpath("blastdb")
        cmd = f"{self.software['makeblastdb']} -in {self.f_fasta} -out {blastdb} -dbtype nucl"
        logger.debug(cmd)
        try:
            subprocess.run(cmd, shell=True, check=True)
        except Exception as e:
            logger.error(e)
            sys.exit(1)

    def _get_primer_pair(self, f_primer):
        """
        获取配对的引物

        :param f_primer: The primer fasta file
        """
        res = []
        for record in SeqIO.parse(f_primer, 'fasta'):
            res.append(record.name)
        paired_res = [(res[i], res[i + 1]) for i in range(0, len(res), 2)]
        return paired_res

    def _parse_blast_result(self, f_blast):
        """
        解析引物的blast结果文件

        :param f_blast: The blast file
        """
        res = {}
        with open(f_blast, 'r') as IN:
            for line in IN:
                if not line.startswith("#"):
                    arr = line.strip().split("\t")
                    res.setdefault(arr[8], [])
                    res[arr[8]].append(arr)
        return res

    def pcr_blast(self, f_primer=None):
        """
        通过blast获取引物扩增的序列
        """
        self.f_primer = Path(f_primer).absolute() if f_primer else self.f_primer

        # blast
        self.d_blast = self.d_out.joinpath("blast")
        self.d_blast.mkdir(exist_ok=True, parents=True)
        ## run blast
        logger.info(f"Use {self.f_primer} blast against {self.f_fasta}")
        f_blast = self.d_blast.joinpath("blast.out")
        cmd = f"{self.software['blastn']} -task blastn-short -db {self.d_prepare.joinpath('blastdb')} -query {self.f_primer} -ungapped -word_size 4 -outfmt '7 sacc sstrand sstart send qlen length nident mismatch qacc qstart qend slen stitle' -out {f_blast}"
        logger.debug(cmd)
        try:
            subprocess.run(cmd, shell=True, check=True)
        except Exception as e:
            logger.error(e)
            sys.exit(1)

        # Get seq info
        info_seq = {}
        for record in SeqIO.parse(self.f_fasta, "fasta"):
            info_seq[record.id] = record.seq

        # Parse the result
        logger.info(f"Parse the blast result")
        f_product = self.d_blast.joinpath("pcr-products.fna")
        handle_product = open(f_product, 'w')
        info_blast = self._parse_blast_result(f_blast)
        info_primer = self._get_primer_pair(self.f_primer)
        for f, r in info_primer:
            info_f = info_blast[f] if f in info_blast else None
            info_r = info_blast[r] if r in info_blast else None
            if info_f and info_r:
                for i in info_f:
                    forward_subject = i[0]
                    forward_strand = i[1]
                    forward_length = int(i[4])
                    forward_ident = int(i[6])
                    if forward_ident + self.mismatches >= forward_length:
                        for j in info_r:
                            reverse_subject = j[0]
                            reverse_strand = j[1]
                            if forward_subject == reverse_subject and forward_strand != reverse_strand:
                                reverse_length = int(j[4])
                                reverse_ident = int(j[6])
                                if reverse_ident + self.mismatches >= reverse_length:
                                    start = int(i[2]) if forward_strand == "plus" else int(j[2])
                                    end = int(j[2]) if forward_strand == "plus" else int(i[2])
                                    len_target = end - start + 1
                                    if self.minlen <= len_target <= self.maxlen:
                                        name = f"{i[8]}|{j[8]}"
                                        seq = info_seq[forward_subject][start - 1:end]
                                        if forward_strand == "minus":
                                            seq = seq.reverse_complement()
                                        handle_product.write(f">{name}\n{seq}\n")
        handle_product.close()

    def get_repeat(self):
        """
        获取重复区域
        """
        logger.info("Extract repeat region")
        f_repeat = self.d_blast.joinpath('repeat-regions.dna')
        cmd = f"""set -e
sed -n '2~2p' {self.d_blast}/pcr-products.fna | while read PRODUCT; do
    REGION="$(echo "$PRODUCT" | grep -Eo '[AG]CA[AC]CAAAA.+..................TA[CT]ATGTCGT' | tail -c +10 | tr -d '\\n' | head -c -28)"
    [ -z "$REGION" ] || echo "$REGION"
done > {f_repeat}
"""
        logger.debug(cmd)
        try:
            subprocess.run(cmd, shell=True, check=True)
        except Exception as e:
            logger.error(e)
            sys.exit(1)

        num_repeat = 0
        with open(f_repeat, 'r') as IN:
            for line in IN:
                if line != "":
                    num_repeat += 1
        if num_repeat == 0:
            logger.error(f"No repeat regions find")
            sys.exit(1)
        if num_repeat > 1:
            logger.warning(f"More than one repeat regions find")

    def type_repeat(self, sample="Query", f_out=None):
        """
        获取分型
        """
        logger.info(f"Type the Sequence")
        self.f_out = Path(f_out).absolute() if f_out else self.f_out
        f_repeat = self.d_blast.joinpath('repeat-regions.dna')
        cmd = f"""set -e
cat {f_repeat} | while read REGION; do
    FOUND="$(grep -Fw "$REGION" {self.db_spa} || grep -Fw "${{REGION%*?}}" {self.db_spa})" || true
    if [ -n "$FOUND" ]; then
        TYPE="$(echo "$FOUND" | cut -d: -f1)"
        FEATURE=$(grep -Fw "$TYPE" {self.type_spa} | cut -d, -f2 || "NA")
        echo -e "{sample}\\t$TYPE\\t$FEATURE"
    fi
done | sort -u > {self.f_out}
sed -i 1i"Sample\\tSpaType\\tRepeatFeature" {self.f_out}
"""
        logger.debug(cmd)
        try:
            subprocess.run(cmd, shell=True, check=True)
        except Exception as e:
            logger.error(e)
            sys.exit(1)

    def run(self):
        """
        Run the pipe
        """
        logger.info(f"Start to run the pipe")
        self.prepare()
        self.pcr_blast()
        self.get_repeat()
        self.type_repeat("Query")

    def to_upload(self, d_upload):
        """
        copy the result file to dir upload
        """
        cmd = f"""set -e
cp -rf {self.f_out} {d_upload}/Spa.tsv
"""
        subprocess.run(cmd, shell=True, check=True)

    @property
    def result(self):
        """
        获取Spa分型结果
        """
        f_result = Path(self.f_out).absolute()
        if f_result.exists():
            with open(f_result, 'r') as IN:
                next(IN)
                for line in IN:
                    arr = line.strip().split("\t")
                    return arr[1]
        else:
            return None


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-n", "--name",
              required=True,
              help="The sample name for your input")
@click.option("-f", "--fasta",
              required=True,
              type=click.Path(),
              help="The input fasta file")
@click.option("--primer",
              required=True,
              type=click.Path(),
              help="The input primer file")
@click.option("--spacsv",
              required=True,
              type=click.Path(),
              help="The spa_db.csv file")
@click.option("--spatype",
              required=True,
              type=click.Path(),
              help="The spatypes.txt file")
@click.option("--mismatches",
              default=2,
              show_default=True,
              type=int,
              help="The max mismatch per primer")
@click.option("--minlen",
              default=100,
              show_default=True,
              type=int,
              help="The minimum required selection length")
@click.option("--maxlen",
              default=5000,
              show_default=True,
              type=int,
              help="The maximum allowed selection length")
@click.option("-o", "--out",
              required=True,
              type=click.Path(),
              help="The out put dir")
def cli(name, fasta, primer, spacsv, spatype, mismatches, minlen, maxlen, out):
    """
    金黄色葡萄球菌Spa型预测
    """
    f_fasta = Path(fasta).absolute()
    f_primer = Path(primer).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)
    f_out = d_out.joinpath(f"{name}.tsv")

    job = Spa(f_fasta, f_primer, spacsv, spatype, d_out, mismatches=mismatches, minlen=minlen, maxlen=maxlen)
    job.prepare()
    job.pcr_blast()
    job.get_repeat()
    job.type_repeat(name, f_out)


if __name__ == "__main__":
    cli()
