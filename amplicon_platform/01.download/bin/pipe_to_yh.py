#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/7/26 13:57
# @Last Modified by:   Ming
# @Last Modified time: 2023/7/26 13:57
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
def Amplicon_Seq_add(row, ud=20):
    chromosome = row["Chrom"]
    seq = info_seq[chromosome]
    start = max(row["Internal_ForwardStart(0base)"] - ud, 0)
    end = min(row["Internal_ReverseEnd(1base)"] + ud, len(seq))
    return str(seq[start:end])


def Amplicon_Seq(row):
    chromosome = row["Chrom"]
    seq = info_seq[chromosome]
    start = row["Internal_ForwardStart(0base)"]
    end = row["Internal_ReverseEnd(1base)"]
    return str(seq[start:end])


def in_primer_specificity(row):
    return f"{row['InternaleF(%)_Specificity']}|{row['InternalR(%)_Specificity']}"


def Amplicon_Seq_add_before(row, up=20):
    chromosome = row["Chrom"]
    seq = info_seq[chromosome]
    end = row["External_ForwardStart(0base)"]
    start = max(end - up, 0)
    return str(seq[start:end])


def External_Amplicon_Seq(row):
    chromosome = row["Chrom"]
    seq = info_seq[chromosome]
    f_start = row["External_ForwardStart(0base)"]
    f_end = row["External_ForwardEnd(1base)"]
    r_start = row["External_ReverseStart(0base)"]
    r_end = row["External_ReverseEnd(1base)"]
    return f"[{seq[f_start:f_end]}]{seq[f_end:r_start]}[{seq[r_start:r_end]}]"


def Amplicon_Seq_add_after(row):
    chromosome = row["Chrom"]
    seq = info_seq[chromosome]
    start = row["External_ReverseEnd(1base)"]
    return seq[start:start + 20]


def out_primer_specificity(row):
    return f"{row['ExternalF(%)_Specificity']}|{row['ExternalR(%)_Specificity']}"


def Block_rev(row):
    chromosome = row["Chrom"]
    seq = info_seq[chromosome]
    f_end = row["External_ForwardEnd(1base)"]
    r_start = row["External_ReverseStart(0base)"]
    return str(seq[f_end:r_start])[::-1]


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-i', '--fin',
              required=True,
              type=click.Path(),
              help="The primer design result from the final pipe")
@click.option('-r', '--ref',
              required=True,
              type=click.Path(),
              help="The reference sequence")
@click.option('-s', '--species',
              required=True,
              type=click.Path(),
              help="The species name")
@click.option('-o', '--fout',
              required=True,
              type=click.Path(),
              help="The final result offer to zhangkai")
def main(fin, ref, species, fout):
    """
    将评估后的最终结果转化为叶辉之前给张开的格式
    """
    f_in = Path(fin).absolute()
    f_ref = Path(ref).absolute()
    f_out = Path(fout).absolute()

    logger.info(f"Read the reference {f_ref}")
    global info_seq
    info_seq = {}
    for record in SeqIO.parse(f_ref, "fasta"):
        name = record.name
        seq = record.seq
        info_seq[name] = seq

    logger.info(f"Read {f_in}")
    df = pd.read_csv(f_in, sep="\t")
    df["species"] = species
    df["ID"] = None
    df["Left_Seq"] = df["Internal_ForwardSeq"]
    df["Right_Seq"] = df["Internal_ReverseSeq"]
    df["Probe_Seq"] = df["Internal_ProbeSeq"]
    df["Amplicon_Name"] = None
    df["Chrom"] = df["Internal_Chromosome"]
    df["Amp_Start"] = df["Internal_ForwardStart(0base)"] + 1
    df["Amp_End"] = df["Internal_ReverseEnd(1base)"]
    df["Amp_Length"] = df["Internal_ReverseEnd(1base)"] - df["Internal_ForwardStart(0base)"]
    df["Left_Tm"] = df["Internal_ForwardTM"]
    df["Right_Tm"] = df["Internal_ReverseTM"]
    df["Probe_Tm"] = df["Internal_ProbeTM"]
    df["Left_Length"] = df["Internal_ForwardLength"]
    df["Right_Length"] = df["Internal_ReverseLength"]
    df["Probe_Length"] = df["Internal_ProbeLength"]
    df["Left_GC"] = df["Internal_ForwardGC"]
    df["Right_GC"] = df["Internal_ReverseGC"]
    df["Probe_GC"] = df["Internal_ProbeGC"]
    df["Taxid"] = df["Species"]
    df["Amplicon_Seq_add"] = df.apply(Amplicon_Seq_add, axis=1)
    df["Amplicon_Seq"] = df.apply(Amplicon_Seq, axis=1)
    df["in_primer_pair_specificity"] = df["InternalPair(%)_Specificity"]
    df["in_primer_specificity"] = df.apply(in_primer_specificity, axis=1)
    df["in_primer_inclusiveness"] = df["InternalPair(%)_Inclusive"]
    df["probe-specificity"] = df["Probe(%)_Specificity"]
    df["probe_inclusiveness"] = df["Probe(%)_Inclusive"]
    df["Internal_USE"] = None
    df["External_ID"] = None
    df["External_Left_Seq"] = df["External_ForwardSeq"]
    df["External_Right_Seq"] = df["External_ReverseSeq"]
    df["External_Amplicon_Name"] = None
    df["External_Chrom"] = df["External_Chromosome"]
    df["External_Amp_Start"] = df["External_ForwardStart(0base)"] + 1
    df["External_Amp_End"] = df["External_ReverseEnd(1base)"]
    df["Amp_Length"] = df["External_ReverseEnd(1base)"] - df["External_ForwardStart(0base)"]
    df["Amplified_Block_Start"] = df["External_ForwardEnd(1base)"] + 1
    df["Amplified_Block_End"] = df["External_ReverseStart(0base)"]
    df["External_Left_Tm"] = df["External_ForwardTM"]
    df["External_Right_Tm"] = df["External_ReverseTM"]
    df["External_Left_Length"] = df["External_ForwardLength"]
    df["External_Right_Length"] = df["External_ReverseLength"]
    df["External_Left_GC"] = df["External_ForwardGC"]
    df["External_Right_GC"] = df["External_ReverseGC"]
    df["Desired_Multiplex"] = None
    df["Designed_Multiplex"] = None
    df["External_Taxid"] = df["Species"]
    df["Amplicon_Seq_add_before"] = df.apply(Amplicon_Seq_add_before, axis=1)
    df["External_Amplicon_Seq"] = df.apply(External_Amplicon_Seq, axis=1)
    df["Amplicon_Seq_add_after"] = df.apply(Amplicon_Seq_add_after, axis=1)
    df["out_primer_pair_specificity"] = df["ExternalPair(%)_Specificity"]
    df["out_primer_specificity"] = df.apply(out_primer_specificity, axis=1)
    df["out_primer_inclusiveness"] = df["ExternalPair(%)_Inclusive"]
    df["Use"] = None
    df["Block_rev"] = df.apply(Block_rev, axis=1)
    df["level"] = None

    logger.info(f"Output to {f_out}")
    columns = ["species", "ID", "Left_Seq", "Right_Seq", "Probe_Seq", "Amplicon_Name", "Chrom", "Amp_Start", "Amp_End",
               "Amp_Length", "Left_Tm", "Right_Tm", "Probe_Tm", "Left_Length", "Right_Length", "Probe_Length",
               "Left_GC", "Right_GC", "Probe_GC", "Taxid", "Amplicon_Seq_add", "Amplicon_Seq",
               "in_primer_pair_specificity", "in_primer_specificity", "in_primer_inclusiveness", "probe-specificity",
               "probe_inclusiveness", "Internal_USE", "External_ID", "External_Left_Seq", "External_Right_Seq",
               "External_Amplicon_Name", "External_Chrom", "External_Amp_Start", "External_Amp_End", "Amp_Length",
               "Amplified_Block_Start", "Amplified_Block_End", "External_Left_Tm", "External_Right_Tm",
               "External_Left_Length", "External_Right_Length", "External_Left_GC", "External_Right_GC",
               "Desired_Multiplex", "Designed_Multiplex", "External_Taxid", "Amplicon_Seq_add_before",
               "External_Amplicon_Seq", "Amplicon_Seq_add_after", "out_primer_pair_specificity",
               "out_primer_specificity", "out_primer_inclusiveness", "Use", "Block_rev", "level"]
    df.to_csv(f_out, sep="\t", columns=columns, index=False, na_rep="NA")


if __name__ == "__main__":
    main()
