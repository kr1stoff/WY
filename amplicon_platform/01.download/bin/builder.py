#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/8/25 10:48
# @Last Modified by:   Ming
# @Last Modified time: 2023/8/25 10:48
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
def get_bac_ref(f_metadata):
    """
    获取细菌的参考基因组
    """
    with open(f_metadata, 'r') as IN:
        for line in IN:
            arr = line.strip().split('\t')
            if arr[4] == "representative genome":
                return arr[0]


def index(d_in: Path, bwa, makeblastdb):
    """
    为ref.fna,refseq.fna,all.fna构建索引

    :param d_in: The input dir contain ref.fna refseq.fna all.fna
    """
    f_ref = d_in.joinpath("ref.fna").absolute()
    f_refseq = d_in.joinpath("refseq.fna").absolute()
    f_refseq_base = d_in.joinpath("refseq")
    f_all = d_in.joinpath("all.fna").absolute()
    f_all_base = d_in.joinpath("all")

    cmd = f"""set -e
{bwa} index {f_ref} >>{f_log} 2>&1
{makeblastdb} -in {f_refseq} -out {f_refseq_base} -dbtype nucl -input_type fasta -title refseq >>{f_log} 2>&1
{makeblastdb} -in {f_all} -out {f_all_base} -dbtype nucl -input_type fasta -title all >>{f_log} 2>&1
"""
    try:
        subprocess.run(cmd, shell=True)
    except:
        sys.exit(logger.error(f"Build index error"))

#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('--taxonid',
              required=True,
              type=click.Path(),
              help="The taxonid's sep by ',' to download")
@click.option('--ptype',
              required=True,
              type=click.Choice(['Bacteria', "Viruses", "Fungi"]),
              help="The pathogen type")
@click.option('--db',
              required=True,
              type=click.Path(),
              help="The database dir")
@click.option('--name',
              required=True,
              type=click.Path(),
              help="The dir name for species database to store")
@click.option('--ref',
              type=click.Path(),
              help="The seq name used for reference")
# @click.option('--latin',
#               type=click.Path(),
#               help="The latin name")
@click.option('--makeblastdb',
              default="makeblastdb",
              type=click.Path(),
              help="The makeblastdb path")
@click.option('--bwa',
              default="bwa",
              type=click.Path(),
              help="The bwa path")
def main(taxonid, db, name, ref, makeblastdb, bwa):
    """
    构建扩增子数据库

    TODO: 根据ptype判断如何建库
    """
    d_database = Path(db).absolute()
    if not d_database.exists():
        sys.exit(logger.error(f"{d_database} not exists!"))

    d_out = d_database.joinpath(name)
    if not d_out.exists():
        sys.exit(logger.error(f"{d_out} not exists!"))
    global f_log
    f_log = d_out.joinpath("build.log")

    d_seq = d_out.joinpath("seq")
    if not d_seq.exists():
        sys.exit(logger.error(f"{d_seq} not exists! Please download the seq first"))

    # ref.fna
    f_ref = d_out.joinpath("ref.fna")
    logger.info(f"Get the {f_ref}")
    if ref:
        # TODO: 病毒的时候应该给什么名字

        cmd = f"zcat {d_seq}/{ref}* > {f_ref}"
    else:
        f_metadata = d_out.joinpath("metadata.tsv")
        ref = get_bac_ref(f_metadata)
        cmd = f"zcat {d_seq}/{ref}* > {f_ref}"
    try:
        subprocess.run(cmd, shell=True)
    except:
        sys.exit(logger.error(f"Can't find the {ref} in your {d_seq}"))

    # refseq.fna
    f_refseq = d_out.joinpath("refseq.fna")
    logger.info(f"Generate the {f_refseq}")
    cmd = f"for i in {d_seq}/GCF*.fna.gz;do zcat $i;done > {f_refseq}"
    try:
        subprocess.run(cmd, shell=True)
    except:
        sys.exit(logger.error(f"Generate {f_refseq} error"))

    # all.fna
    # TODO: 添加质控
    f_all = d_out.joinpath("all.fna")
    logger.info(f"Generate the {f_all}")
    cmd = f"for i in {d_seq}/*.fna.gz;do zcat $i;done > {f_all}"
    try:
        subprocess.run(cmd, shell=True)
    except:
        sys.exit(logger.error(f"Generate {f_all} error"))

    # 建库
    logger.info(f"Build the index")
    index(d_out, bwa, makeblastdb)

    # assembly_chromosome.tsv
    d_bin = Path(__file__).absolute().parent
    f_assembly = d_out.joinpath("assembly_chromosome.tsv")
    logger.info(f"Generate the {f_assembly}")
    cmd = f"python {d_bin}/assembly_chromosome.py --dirin {d_seq} -t {taxonid} -o {f_assembly}"
    try:
        subprocess.run(cmd, shell=True)
    except:
        sys.exit(logger.error(f"Generate {f_assembly} error"))


if __name__ == "__main__":
    main()
