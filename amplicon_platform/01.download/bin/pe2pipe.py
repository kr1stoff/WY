#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/2/15 9:45
# @Last Modified by:   Ming
# @Last Modified time: 2023/2/15 9:45
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
@click.option('--pe',
              required=True,
              type=click.Path(),
              help="The PrimerExpress3.1 result")
@click.option("--sequence",
              required=True,
              type=click.Path(),
              help="The region sequence")
@click.option("--region",
              required=True,
              type=click.Path(),
              help="The region bed file indicate it's position in the genome/chromosome")
@click.option("--name",
              required=True,
              help="The region name")
@click.option('--prefix',
              required=True,
              type=click.Path(),
              help="The output prefix")
def main(pe, sequence, region, name, prefix):
    """
    将PrimerExpress3的结果转为流程内引物预测结果
    """
    f_pe = Path(pe).absolute()
    f_sequence = Path(sequence).absolute()
    f_region = Path(region).absolute()

    d_out = Path(prefix).absolute().parent
    d_out.mkdir(exist_ok=True)
    prefix = Path(prefix).absolute()

    logger.info(f"Parse the sequence info")
    info_seq = {}
    for record in SeqIO.parse(f_sequence, "fasta"):
        info_seq[record.name] = record.seq

    logger.info(f"Parse the region info {f_region}")
    info_region = {}
    with open(f_region, 'r') as IN:
        for line in IN:
            arr = line.strip().split("\t")
            info_region[arr[3]] = arr

    logger.info(f"Parse the PrimerExpress3.1 result")
    region_name = name
    chromosome_name = info_region[region_name][0]
    region_start = int(info_region[region_name][1])
    region_seq = info_seq[region_name]
    df = pd.read_csv(f_pe, sep="\t")

    # Start转为0base
    df["Fwd Start"] = df["Fwd Start"] - 1
    df["Rev Start"], df["Rev Stop"] = df["Rev Stop"], df["Rev Start"]
    df["Rev Start"] = df["Rev Start"] - 1
    df["Probe Start"] = df["Probe Start"] - 1
    # 获取包含大小写的序列
    df["ForwardSeq"] = df.apply(lambda row: region_seq[row["Fwd Start"]:row["Fwd Stop"]],
                                axis=1)
    df["ReverseSeq"] = df.apply(lambda row: region_seq[row["Rev Start"]:row["Rev Stop"]].reverse_complement(),
                                axis=1)
    df["ProbeSeq"] = df.apply(lambda row: region_seq[row["Probe Start"]:row["Probe Stop"]],
                              axis=1)
    df.rename(columns={'#': "Index",
                       # "Fwd Start": "ForwardStart(0base)",
                       # "Fwd Stop": "ReverseEnd(1base)",
                       "Fwd Length": "ForwardLength",
                       "Fwd Tm": "ForwardTM",
                       "Fwd %GC": "ForwardGC",
                       # "Fwd Seq": "ForwardSeq",
                       # "Rev Start": "ReverseStart(0base)",
                       # "Rev Stop": "ReverseEnd(1base)",
                       "Rev Length": "ReverseLength",
                       "Rev Tm": "ReverseTM",
                       "Rev %GC": "ReverseGC",
                       # "Rev Seq": "ReverseSeq",
                       # "Probe Start": "ProbeStart(0base)",
                       # "Probe Stop": "ProbeEnd(1base)",
                       "Probe Length": "ProbeLength",
                       "Probe Tm": "ProbeTM",
                       "Probe %GC": "ProbeGC",
                       # "Probe Seq": "ProbeSeq"
                       },
              inplace=True)
    df["Chromosome"] = chromosome_name
    # 转为对应的参考基因组上的坐标
    df["ForwardStart(0base)"] = df["Fwd Start"] + region_start
    df["ForwardEnd(1base)"] = df["Fwd Stop"] + region_start
    df["ReverseStart(0base)"] = df["Rev Start"] + region_start
    df["ReverseEnd(1base)"] = df["Rev Stop"] + region_start
    df["ProbeStart(0base)"] = df["Probe Start"] + region_start
    df["ProbeEnd(1base)"] = df["Probe Stop"] + region_start

    f_info = f"{prefix}.xls"
    logger.info(f"Output the info to {f_info}")
    df.to_csv(f_info, sep="\t", index=False,
              columns=["Index", "Chromosome",
                       "ForwardStart(0base)", "ForwardEnd(1base)", "ForwardLength", "ForwardTM", "ForwardGC",
                       "ForwardSeq",
                       "ReverseStart(0base)", "ReverseEnd(1base)", "ReverseLength", "ReverseTM", "ReverseGC",
                       "ReverseSeq",
                       "ProbeStart(0base)", "ProbeEnd(1base)", "ProbeLength", "ProbeTM", "ProbeGC", "ProbeSeq"])
    f_bed = f"{prefix}.bed"
    logger.info(f"Output the bed info to {f_bed}")
    with open(f_bed, 'w') as OUT:
        for i in range(df.shape[0]):
            print(*[df.iloc[i]["Chromosome"],
                    df.iloc[i]["ForwardStart(0base)"],
                    df.iloc[i]["ForwardEnd(1base)"],
                    f"{df.iloc[i]['Chromosome']}_{df.iloc[i]['Index']}-forward",
                    0, "+"],
                  sep="\t", file=OUT)
            print(*[df.iloc[i]["Chromosome"],
                    df.iloc[i]["ReverseStart(0base)"],
                    df.iloc[i]["ReverseEnd(1base)"],
                    f"{df.iloc[i]['Chromosome']}_{df.iloc[i]['Index']}-reverse",
                    0, "-"],
                  sep="\t", file=OUT)
            print(*[df.iloc[i]["Chromosome"],
                    df.iloc[i]["ProbeStart(0base)"],
                    df.iloc[i]["ProbeEnd(1base)"],
                    f"{df.iloc[i]['Chromosome']}_{df.iloc[i]['Index']}-probe",
                    0, "+"],
                  sep="\t", file=OUT)

    f_fasta = f"{prefix}.fasta"
    logger.info(f"Output the seq info to {f_fasta}")
    with open(f_fasta, 'w') as OUT:
        for i in range(df.shape[0]):
            print(f">{df.iloc[i]['Chromosome']}_{df.iloc[i]['Index']}-forward\n{df.iloc[i]['ForwardSeq']}",
                  file=OUT)
            print(f">{df.iloc[i]['Chromosome']}_{df.iloc[i]['Index']}-reverse\n{df.iloc[i]['ReverseSeq']}",
                  file=OUT)
            print(f">{df.iloc[i]['Chromosome']}_{df.iloc[i]['Index']}-probe\n{df.iloc[i]['ProbeSeq']}",
                  file=OUT)


if __name__ == "__main__":
    main()
