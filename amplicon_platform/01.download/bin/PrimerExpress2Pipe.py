#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/13 16:29
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/13 16:29
import logging
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-pe', '--primerexpress',
              required=True,
              type=click.Path(),
              help="The PrimerExpress output table file")
@click.option('-fa', '--fasta',
              required=True,
              type=click.Path(),
              help="The soft mask fasta name")
@click.option('-n', '--name',
              required=True,
              type=click.Path(),
              help="The Chromosome name")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
def main(primerexpress, fasta, name, out):
    """
    将PrimerExpress3.1的输出结果整理成流程所需的结果
    """
    f_in = Path(primerexpress).absolute()
    f_fasta = Path(fasta).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)

    logger.info(f"Read the PrimerExpress3.1 result table file")
    df = pd.read_table(f_in, sep="\t")
    df.rename(columns={"#": "Index",
                       "Fwd Start": "ForwardStart(0base)",
                       "Fwd Stop": "ForwardEnd(1base)",
                       "Fwd Length": "ForwardLength",
                       "Fwd Tm": "ForwardTM",
                       "Fwd %GC": "ForwardGC",
                       "Fwd Seq": "ForwardSeq",
                       "Rev Start": "ReverseEnd(1base)",
                       "Rev Stop": "ReverseStart(0base)",
                       "Rev Length": "ReverseLength",
                       "Rev Tm": "ReverseTM",
                       "Rev %GC": "ReverseGC",
                       "Rev Seq": "ReverseSeq",
                       "Probe Start": "ProbeStart(0base)",
                       "Probe Stop": "ProbeEnd(1base)",
                       "Probe Length": "ProbeLength",
                       "Probe Tm": "ProbeTM",
                       "Probe %GC": "ProbeGC",
                       "Probe Seq": "ProbeSeq"}, inplace=True)
    df["Chromosome"] = name
    df["ForwardStart(0base)"] = df["ForwardStart(0base)"] - 1
    df["ReverseStart(0base)"] = df["ReverseStart(0base)"] - 1
    df["ProbeStart(0base)"] = df["ProbeStart(0base)"] - 1

    logger.info(f"Get the Chromosome {name} seq from file {f_fasta}")
    target_seq = None
    for record in SeqIO.parse(f_fasta, "fasta"):
        if record.name == name:
            target_seq = record.seq

    logger.info(f"Fetch the primer and probe sequence")
    df["ForwardSeq"] = df.apply(lambda x: str(target_seq[x["ForwardStart(0base)"]:x["ForwardEnd(1base)"]]), axis=1)
    df["ReverseSeq"] = df.apply(
        lambda x: str(target_seq[x["ReverseStart(0base)"]:x["ReverseEnd(1base)"]].reverse_complement()), axis=1)
    df["ProbeSeq"] = df.apply(lambda x: str(target_seq[x["ProbeStart(0base)"]:x["ProbeEnd(1base)"]]), axis=1)

    logger.info(f"Output the result to dir {d_out}")
    f_table = d_out.joinpath("primer.xls")
    columns = ["Index", "Chromosome", "ForwardStart(0base)", "ForwardEnd(1base)", "ForwardLength", "ForwardTM",
               "ForwardGC", "ForwardSeq", "ReverseStart(0base)", "ReverseEnd(1base)", "ReverseLength", "ReverseTM",
               "ReverseGC", "ReverseSeq", "ProbeStart(0base)", "ProbeEnd(1base)", "ProbeLength", "ProbeTM", "ProbeGC",
               "ProbeSeq"]
    df.to_csv(f_table, sep="\t", columns=columns, index=False)
    f_fasta = d_out.joinpath("primer.fasta")
    f_bed = d_out.joinpath("primer.bed")
    with open(f_fasta, 'w') as OUT1, open(f_bed, 'w') as OUT2:
        for row in df.loc[:, columns].iterrows():
            arr = row[1]
            print(f">{arr[1]}_{arr[0]}-forward\n{arr[7]}", file=OUT1)
            print(f">{arr[1]}_{arr[0]}-reverse\n{arr[13]}", file=OUT1)
            print(f">{arr[1]}_{arr[0]}-probe\n{arr[19]}", file=OUT1)
            print(*[arr[1], arr[2], arr[3], f"{arr[1]}_{arr[0]}-forward", 0, "+"], sep="\t", file=OUT2)
            print(*[arr[1], arr[8], arr[9], f"{arr[1]}_{arr[0]}-reverse", 0, "-"], sep="\t", file=OUT2)
            print(*[arr[1], arr[14], arr[15], f"{arr[1]}_{arr[0]}-probe", 0, "+"], sep="\t", file=OUT2)


if __name__ == "__main__":
    main()
