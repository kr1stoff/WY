#!/usr/bin/python3
# @auther: stone
# @date: 2023-08-15
import os
import click
import re
import math
import subprocess

# == Option =====================================================================


def format_number(n, width):
    return str(n).zfill(width)


def num_set(start, end):  # 数字自动扩宽，从01到11
    max_width = len(str(end))
    string_list = []
    for n in range(start, end + 1):
        formatted_number = format_number(n, max_width)
        string_list.append(formatted_number)
    return string_list


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-f', '--fasta', required=True, type=click.Path(), help="输入fa")
@click.option('-o', '--out', default='blastn.out', type=click.Path(), help="输出结果")
@click.option('-n', '--nt_database', default='/sdbb/bioinfor/yehui/nt/nt.fa', type=click.Path(), help="指定blast比对库")
@click.option('-t/', '--threads', default='50', type=str, help="blast线程数")
@click.option('-q/', '--qsub_num', default='3', type=str, help="fa过多时,多少个blast同时跑")
@click.option('-c/', '--cuts', default='100', type=int, help="切成的")
@click.option('--short/--long', default=True, show_default=True, help="短序列还是长序列")
def main(fasta, out, nt_database, short, threads, qsub_num, cuts):
    Bin = os.path.dirname(os.path.abspath(__file__))
    fasta = os.path.abspath(fasta)
    fasta_name = os.path.basename(fasta)
    out = os.path.abspath(out)
    outdir = os.path.dirname(out)
    nt_database = os.path.abspath(nt_database)
    contig_num = 0
    with open(fasta, 'r') as IN:
        for line in IN:
            if re.match(r'^>', line):
                contig_num += 1
    cutf = math.ceil(contig_num / cuts)
    format_nums = num_set(1, cutf)
    if (short):
        opt = '-task blastn-short -word_size 9 -evalue 800000 -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1'
        outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs staxid sscinames sstrand'"

    else:
        opt = '-evalue 1e-5 '
        outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus staxid sscinames scomnames'"

    if (contig_num >= cuts):
        with open(f'{outdir}/blastn.sh', 'w') as BLAST_SHELL:
            blast_shell = ''
            for i in format_nums:
                blast_shell += f"""export BLASTDB=/sdbb/bioinfor/yehui/nt
blastn -query {outdir}/{fasta_name}.cut/{fasta_name}.{i} -db {nt_database} -num_threads {threads} {opt} -outfmt {outfmt} -out {out}.{i}
"""
            print(blast_shell, file=BLAST_SHELL)
        shell = f"""
perl {Bin}/fastaDeal.pl --cuts {cuts} {fasta} --outdir {outdir}
rm {out}.* -rf
perl {Bin}/qsub_all_local.pl -b 2 -m {qsub_num} -d {outdir}/qsub_blastn {outdir}/blastn.sh
cat {out}.* > {out}
rm {out}.* {outdir}/{fasta_name}.cut/ {outdir}/qsub_blastn -rf
"""
    else:
        shell = f"""export BLASTDB=/sdbb/bioinfor/yehui/nt
blastn -query {fasta} -db {nt_database} -num_threads {threads} {opt} -outfmt {outfmt} -out {out}
"""
    subprocess.call(shell, shell=True)
    # print(shell)


if __name__ == "__main__":
    main()
