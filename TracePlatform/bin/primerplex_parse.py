#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/2/13 10:05
# @Last Modified by:   Ming
# @Last Modified time: 2023/2/13 10:05
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
@click.option('--primerplex',
              required=True,
              type=click.Path(),
              help="The primerplex info file")
@click.option('--ref',
              required=True,
              type=click.Path(),
              help="The ref sequence")
@click.option('-p', '--prefix',
              required=True,
              type=click.Path(),
              help="The out put prefix")
def main(primerplex, ref, prefix):
    """
    处理PrimerPlex的结果
    """
    f_primerplex = Path(primerplex).absolute()
    f_ref = Path(ref).absolute()
    pfx = Path(prefix).name
    d_out = Path(prefix).absolute().parent
    d_out.mkdir(exist_ok=True)

    logger.info(f"Parse the ref info")
    for record in SeqIO.parse(f_ref, "fasta"):
        info_seq = record

    logger.info(f"Read the PrimerPlex result file {f_primerplex}")
    df = pd.read_excel(f_primerplex)
    df.rename(columns={'#': "Index",
                       "Chrom": "Chromosome",
                       "Amplicon_Name": "InternalRegion",
                       "Amplicon_Start": "ForwardStart(0base)",
                       "Amplicon_End": "ReverseEnd(1base)",
                       "Left_Primer_Tm": "ForwardTM",
                       "Right_Primer_Tm": "ReverseTM",
                       "Left_Primer_Length": "ForwardLength",
                       "Right_Primer_Length": "ReverseLength",
                       "Left_GC": "ForwardGC",
                       "Right_GC": "ReverseGC",
                       "Left_Primer_Length": "ForwardLength",
                       "Right_Primer_Length": "ReverseLength"},
              inplace=True)
    df["ForwardEnd(1base)"] = df["ForwardStart(0base)"] + df["ForwardLength"]
    df["ReverseEnd(1base)"] = df["ReverseEnd(1base)"] + 1
    df["ReverseStart(0base)"] = df["ReverseEnd(1base)"] - df["ReverseLength"]
    df["ForwardSeq"] = [str(info_seq[start:end].seq) for start, end in
                        zip(df["ForwardStart(0base)"], df["ForwardEnd(1base)"])]
    df["ReverseSeq"] = [str(info_seq[start:end].seq.reverse_complement()) for start, end in
                        zip(df["ReverseStart(0base)"], df["ReverseEnd(1base)"])]

    f_info = d_out.joinpath(f"{pfx}.xls")
    logger.info(f"Output the info to {f_info}")
    df.to_csv(f_info, sep="\t", index=False,
              columns=["Index", "Chromosome", "InternalRegion", "ForwardStart(0base)", "ForwardEnd(1base)",
                       "ForwardLength", "ForwardTM", "ForwardGC", "ForwardSeq", "ReverseStart(0base)",
                       "ReverseEnd(1base)", "ReverseLength", "ReverseTM", "ReverseGC", "ReverseSeq"])

    f_bed = d_out.joinpath(f"{pfx}.bed")
    logger.info(f"Output the bed info to {f_bed}")
    with open(f_bed, 'w') as OUT:
        for i in range(df.shape[0]):
            name = f"{df.iloc[i]['InternalRegion']}"
            print(*[df.iloc[i]['Chromosome'], df.iloc[i]['ForwardStart(0base)'], df.iloc[i]['ForwardEnd(1base)'],
                    f"{name}_external_forward", 0, "+"], sep="\t", file=OUT)
            print(*[df.iloc[i]['Chromosome'], df.iloc[i]['ReverseStart(0base)'], df.iloc[i]['ReverseEnd(1base)'],
                    f"{name}_external_reverse", 0, "-"], sep="\t", file=OUT)

    f_fa = d_out.joinpath(f"{pfx}.fa")
    logger.info(f"Output the fa info to {f_fa}")
    with open(f_fa, 'w') as OUT:
        for i in range(df.shape[0]):
            name = f"{df.iloc[i]['InternalRegion']}"
            print(f">{name}_external_forward\n{df.iloc[i]['ForwardSeq']}", file=OUT)
            print(f">{name}_external_reverse\n{df.iloc[i]['ReverseSeq']}", file=OUT)


if __name__ == "__main__":
    main()
