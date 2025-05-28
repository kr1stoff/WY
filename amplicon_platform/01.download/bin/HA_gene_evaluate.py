#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/6/25 17:09
# @Last Modified by:   Ming
# @Last Modified time: 2023/6/25 17:09
import logging
from pathlib import Path

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def blast(f_fasta, database, f_out, blastn, param=None):
    """
    Get the blast cmd
    """
    param = "-num_threads 60 -word_size 9 -evalue 800000 -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sstrand'" if param is None else param
    cmd = f"{blastn} -task blastn-short -query {f_fasta} -db {database} {param} -out {f_out}"
    return cmd


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option("--primer",
              required=True,
              type=click.Path(),
              help="The primer fasta file")
@click.option("--primer_name",
              required=True,
              type=click.Path(),
              help="The primer name info file")
@click.option("--probe_name",
              required=False,
              type=click.Path(),
              help="The probe name info file")
@click.option("-db", "--database",
              required=True,
              type=click.Path(),
              help="The H gene sequence blast database")
@click.option("--info",
              required=True,
              type=click.Path(),
              help="The info file contain the HA segment number in database")
@click.option("--out",
              required=True,
              type=click.Path(),
              help="The out put dir")
@click.option("--blastn",
              show_default=True,
              default="/home/yehui/software/bin/blastn",
              type=click.Path(),
              help="The blastn path")
@click.option("--perl",
              default="/usr/bin/perl",
              show_default=True,
              type=click.Path(),
              help="The perl path")
@click.option("--python",
              show_default=True,
              default="/sdbb/bioinfor/renchaobo/Software/miniconda3/bin/python",
              type=click.Path(),
              help="The python path")
def main(primer, primer_name, probe_name, database, info, out, blastn, perl, python):
    """
    评估甲流H基因不同分型的特异性和包容性
    """
    f_primer = Path(primer).absolute()
    f_primer_name = Path(primer_name).absolute()
    f_probe_name = Path(probe_name).absolute() if probe_name else None
    database = Path(database).absolute()
    f_info = Path(info).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)
    d_bin = Path(__file__).absolute().parent

    logger.info(f"Start to Generate the command")
    cmd = ["set -e"]
    # Blast
    d_blast = d_out.joinpath("Blast")
    d_blast.mkdir(exist_ok=True, parents=True)
    f_blast_out = d_blast.joinpath("blast.tsv")
    cmd.append(blast(f_primer, database, f_blast_out, blastn))
    # 引物评估
    d_evaluate = d_out.joinpath("Evaluate")
    d_evaluate.mkdir(exist_ok=True, parents=True)
    ## 过滤
    f_blast_out_filter = d_evaluate.joinpath("blast.filter.tsv")
    cmd.append(f"{perl} {d_bin}/filter_blast.pl {f_blast_out} > {f_blast_out_filter}")
    ## 引物特异性检测
    f_primer_pair_evaluate = d_evaluate.joinpath("primer.stat.xls")
    cmd.append(
        f"{python} {d_bin}/HA_primer_pair_evaluate.py --primer_name {f_primer_name} --blast {f_blast_out_filter}  --info {f_info} --out {f_primer_pair_evaluate}")
    ## 探针特异性检测
    if f_probe_name:
        f_probe_evaluate = d_evaluate.joinpath("probe.stat.xls")
        cmd.append(
            f"{python} {d_bin}/HA_probe_evaluate.py --probe_name {f_probe_name} --blast {f_blast_out_filter} --info {f_info} --out {f_probe_evaluate}")

    f_script = d_out.joinpath("run.sh")
    logger.info(f"Write the command to script: {f_script}")
    with open(f_script, 'w') as OUT:
        print(*cmd, sep="\n", file=OUT)


if __name__ == "__main__":
    main()
