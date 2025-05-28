#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/8/8 16:56
# @Last Modified by:   Ming
# @Last Modified time: 2023/8/8 16:56
import logging
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def Amplicon_Seq(row):
    """

    """
    chromosome = row["Chromosome"]
    start = row["ForwardStart(0base)"]
    end = row["ReverseEnd(1base)"]
    return info_seq[chromosome][start:end]


def Primer_perecnt(row):
    """

    """
    res = row["Inclusive Pair Number(percent %)"].strip().split('(')[1].rstrip(')')
    return res


def Primer_total(row):
    """

    """
    res = row["Inclusive Pair Number(percent %)"].strip().split('(')[0]
    return res


def Genome_num(row):
    """

    """
    percent = float(row["Inclusive Pair Number(percent %)"].strip().split('(')[1].rstrip(')'))
    num = int(row["Inclusive Pair Number(percent %)"].strip().split('(')[0])
    res = round(num / (percent / 100))
    return res


def Pair_percent(row):
    """

    """
    percent = float(row["Specificity Pair Number(percent %)"].strip().split('(')[1].rstrip(')'))
    return percent


def Pair_total(row):
    """

    """
    num = int(row["Specificity Pair Number(percent %)"].strip().split('(')[0])
    return num


def F_percent(row):
    """
    """
    percent = float(row["Specificity F Number(percent %)"].strip().split('(')[1].rstrip(')'))
    return percent


def F_total(row):
    """

    """
    num = int(row["Specificity F Number(percent %)"].strip().split('(')[0])
    return num


def R_percent(row):
    """

    """
    percent = float(row["Specificity R Number(percent %)"].strip().split('(')[1].rstrip(')'))
    return percent


def R_total(row):
    """

    """
    num = int(row["Specificity R Number(percent %)"].strip().split('(')[0])
    return num


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-i', '--fin',
              required=True,
              type=click.Path(),
              help="The input file")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output file")
@click.option('-r', '--ref',
              required=True,
              type=click.Path(),
              help="The reference seq file")
def main(fin, out, ref):
    """
    流程结果转为给张开的格式
    """
    f_in = Path(fin).absolute()
    f_out = Path(out).absolute()
    f_ref = Path(ref).absolute()

    logger.info(f"Parse the reference {f_ref}")
    global info_seq
    info_seq = {}
    for record in SeqIO.parse(f_ref, "fasta"):
        info_seq[record.name] = record.seq

    logger.info(f"Read the {f_in}")
    df = pd.read_csv(f_in, sep="\t")
    df["ID"] = None
    df["Left_Seq"] = df["ForwardSeq"]
    df["Right_Seq"] = df["ReverseSeq"]
    df["Amplicon_Name"] = df["TargetRegion"]
    df["Chrom"] = df["Chromosome"]
    df["Amp_Start"] = df["ForwardStart(0base)"] + 1
    df["Amp_End"] = df["ReverseEnd(1base)"]
    df["Amp_Length"] = df["ReverseEnd(1base)"] - df["ForwardStart(0base)"]
    df["Amplified_Block_Start"] = df["ForwardEnd(1base)"] + 1
    df["Amplified_Block_End"] = df["ReverseStart(0base)"]
    df["Left_Tm"] = df["ForwardTM"]
    df["Right_Tm"] = df["ReverseTM"]
    df["Left_Length"] = df["ForwardLength"]
    df["Right_Length"] = df["ReverseLength"]
    df["Left_GC"] = df["ForwardGC"]
    df["Right_GC"] = df["ReverseGC"]
    df["Desired_Multiplex"] = None
    df["Designed_Multiplex"] = None
    df["Taxid"] = df["Species"]
    df["Amplicon_Seq"] = df.apply(Amplicon_Seq, axis=1)
    # 特异性
    df["Pair_perecnt"] = df.apply(Pair_percent, axis=1) if "Specificity Pair Number(percent %)" in df.columns else None
    df["Pair_total"] = df.apply(Pair_total, axis=1) if "Specificity Pair Number(percent %)" in df.columns else None
    df["F_percent"] = df.apply(F_percent, axis=1) if "Specificity F Number(percent %)" in df.columns else None
    df["F_total"] = df.apply(F_total, axis=1) if "Specificity F Number(percent %)" in df.columns else None
    df["R_percent"] = df.apply(R_percent, axis=1) if "Specificity R Number(percent %)" in df.columns else None
    df["R_total"] = df.apply(R_total, axis=1) if "Specificity R Number(percent %)" in df.columns else None
    # 包容性
    df["Primer_perecnt"] = df.apply(Primer_perecnt, axis=1)
    df["Primer_total"] = df.apply(Primer_total, axis=1)
    df["Genome_num"] = df.apply(Genome_num, axis=1)
    df["Use"] = None

    logger.info(f"Output to {f_out}")
    if "Specificity Pair Number(percent %)" in df.columns:
        columns = ["ID", "Left_Seq", "Right_Seq", "Amplicon_Name", "Chrom", "Amp_Start", "Amp_End", "Amp_Length",
                   "Amplified_Block_Start", "Amplified_Block_End", "Left_Tm", "Right_Tm", "Left_Length", "Right_Length",
                   "Left_GC", "Right_GC", "Desired_Multiplex", "Designed_Multiplex", "Taxid", "Amplicon_Seq",
                   "Pair_perecnt", "Pair_total", "F_percent", "F_total", "R_percent", "R_total",
                   "Primer_perecnt", "Primer_total", "Genome_num", "Use"]
    else:
        columns = ["ID", "Left_Seq", "Right_Seq", "Amplicon_Name", "Chrom", "Amp_Start", "Amp_End", "Amp_Length",
                   "Amplified_Block_Start", "Amplified_Block_End", "Left_Tm", "Right_Tm", "Left_Length", "Right_Length",
                   "Left_GC", "Right_GC", "Desired_Multiplex", "Designed_Multiplex", "Taxid", "Amplicon_Seq",
                   "Primer_perecnt", "Primer_total", "Genome_num", "Use"]
    df.to_csv(f_out, sep="\t", columns=columns, index=False, na_rep="NA")


if __name__ == "__main__":
    main()
