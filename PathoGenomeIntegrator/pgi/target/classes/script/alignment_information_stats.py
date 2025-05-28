#!/usr/bin/env python

import re
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import click

matplotlib.use('Agg')


@click.command()
@click.option('-a', '--accessions', required=True, type=click.Path(exists=True), help='输入 Accession 文件')
@click.option('-A', '--analysis', required=True, type=click.Path(exists=True), help='输入分析结果目录')
@click.option('-o', '--output', required=True, type=click.Path(exists=False), help='输出结果目录')
def main(accessions, analysis, output):

    analysis_dir = Path(analysis)
    output_dir = Path(output)
    df = get_accession_reference_parameter(accessions, analysis_dir, output_dir)
    draw_stats_figure(df, output_dir)


def get_accession_reference_parameter(acces, anal_dir: Path, outdir: Path):
    df = pd.DataFrame()

    with open(acces) as f:
        for line in f:
            acce = line.strip()
            df = concat_new_record(anal_dir.joinpath(
                f'align_fastq_to_integrated/{acce}/{acce}.bam.stats'), acce, 'integrate', df)
            df = concat_new_record(anal_dir.joinpath(
                f'align_fastq_to_reference/{acce}/{acce}.bam.stats'), acce, 'reference', df)

    df.columns = ['accession', 'genome', 'parameter', 'value']
    df.to_csv(outdir.joinpath('accession_reference_parameter.tsv'), index=False, sep='\t')
    return df


def get_mapped_rate(file: str):
    dicbs = {}
    with open(file) as f:
        for line in f:
            ls = re.sub('\t#.*', '', line).strip().split(':\t')
            dicbs[ls[0]] = ls[1]
    read_mapped_rate = int(dicbs['reads mapped']) / int(dicbs['raw total sequences'])
    base_mapped_rate = int(dicbs['bases mapped (cigar)']) / \
        int(dicbs['total length'])  # cigar base mapped (more accurate)
    error_rate = float(dicbs['error rate'])  # mismatches / bases mapped (cigar)
    return read_mapped_rate, base_mapped_rate, error_rate


def concat_new_record(file_bam_stats: str, accession: str, reference: str, datafram: pd.DataFrame):
    read_mapped_rate, base_mapped_rate, error_rate = get_mapped_rate(file_bam_stats)
    datafram = pd.concat([datafram, pd.DataFrame(
        [[accession, reference, 'read_mapped_rate', read_mapped_rate]])], axis=0, ignore_index=True)
    datafram = pd.concat([datafram, pd.DataFrame(
        [[accession, reference, 'base_mapped_rate', base_mapped_rate]])], axis=0, ignore_index=True)
    datafram = pd.concat([datafram, pd.DataFrame(
        [[accession, reference, 'error_rate', error_rate]])], axis=0, ignore_index=True)
    return datafram


def draw_stats_figure(df, output_dir: Path):
    # Mutliple Figures
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(30, 15))
    sns.stripplot(data=df[df['parameter'] == 'read_mapped_rate'],
                  x='genome', y='value', palette='Set2', hue='genome', linewidth=1, legend=False, ax=axes[0, 0])
    sns.stripplot(data=df[df['parameter'] == 'base_mapped_rate'],
                  x='genome', y='value', palette='Set2', hue='genome', linewidth=1, legend=False, ax=axes[0, 1])
    sns.stripplot(data=df[df['parameter'] == 'error_rate'],
                  x='genome', y='value', palette='Set2', hue='genome', linewidth=1, legend=False, ax=axes[0, 2])

    sns.lineplot(data=df[(df['parameter'] == 'read_mapped_rate') & (df['genome'] == 'integrate')],
                 x='accession', y='value', linestyle='-', marker='o', markersize=8,
                 label='Integrate', color='blue', ax=axes[1, 0])
    sns.lineplot(data=df[(df['parameter'] == 'read_mapped_rate') & (df['genome'] == 'reference')],
                 x='accession', y='value', linestyle='--', marker='s', markersize=8,
                 label='Reference', color='green', ax=axes[1, 0])

    sns.lineplot(data=df[(df['parameter'] == 'base_mapped_rate') & (df['genome'] == 'integrate')],
                 x='accession', y='value', linestyle='-', marker='o', markersize=8,
                 label='Integrate', color='blue', ax=axes[1, 1])
    sns.lineplot(data=df[(df['parameter'] == 'base_mapped_rate') & (df['genome'] == 'reference')],
                 x='accession', y='value', linestyle='--', marker='s', markersize=8,
                 label='Reference', color='green', ax=axes[1, 1])

    sns.lineplot(data=df[(df['parameter'] == 'error_rate') & (df['genome'] == 'integrate')],
                 x='accession', y='value', linestyle='-', marker='o', markersize=8,
                 label='Integrate', color='blue', ax=axes[1, 2])
    sns.lineplot(data=df[(df['parameter'] == 'error_rate') & (df['genome'] == 'reference')],
                 x='accession', y='value', linestyle='--', marker='s', markersize=8,
                 label='Reference', color='green', ax=axes[1, 2])

    for ax in axes[0]:
        ax.set_ylabel('Rate')
        ax.set_xlabel('Genome')

    for ax in axes[1]:
        ax.set_ylabel('Rate')
        ax.set_xlabel('Accession')
        ax.set_xticklabels([])

    axes[0, 0].set_title('Reads Mapped Rate Across Genomes')
    axes[0, 1].set_title('Bases Mapped Rate Across Genomes')
    axes[0, 2].set_title('Error Rate Across Genomes')
    axes[1, 0].set_title('Reads Mapped Rate Across Accessions')
    axes[1, 1].set_title('Bases Mapped Rate Across Accessions')
    axes[1, 2].set_title('Error Rate Across Accessions')

    plt.savefig(output_dir.joinpath('multiple_stats_diff.png'), dpi=300)


if __name__ == '__main__':
    main()
