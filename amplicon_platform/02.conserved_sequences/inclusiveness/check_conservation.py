#!/usr/bin/python3
# @auther: zhuzi
# @date: 2022-12-05
from subprocess import PIPE, run, Popen
import sys
import os
import click

p_tools = "/sdbb/bioinfor/zhuzi/tools"

# == Option =====================================================================


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("bed", required=True, type=click.Path())
@click.option('-f', '--fasta', required=True, type=click.Path(), help="原始基因集，用于提取bed中的fa")
@click.option('-db', required=True, type=click.Path(), help="指定blast比对库")
@click.option('-c', '--chrome_file', required=True, type=click.Path(), help="染色体文件")
@click.option('-o', '--out', default='./', type=click.Path(), help="输出目录,默认[./]")
def main(bed, fasta, out, chrome_file, db):
    """评估bed文件基因的保守性，\n[bed format]: gene_name\\tstart\\tend"""
    out = os.path.abspath(out)
    os.makedirs(out, exist_ok=True)
    bed = os.path.abspath(bed)
    fasta = os.path.abspath(fasta)
    chrome_file = os.path.abspath(chrome_file)
    db = os.path.abspath(db)

    cmd = f"""
# fetch fasta
bedtools getfasta -fi {fasta} -bed {bed} -fo {out}/0.bed.fa

# blast
blastn -db {db} -evalue 1e-5 -num_threads 40 -max_target_seqs 2000 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus' -query {out}/0.bed.fa -out {out}/0.blast

# filter --cov 95 --identity 95
python3 {p_tools}/filter_blast.py filter  {out}/0.blast -f 15,3 --cov 95 --identity 95 -o {out}/1.blast.filter

# add-header
sed -i "1iqaccver\\tsaccver\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tqcovs\\tqcovhsp\\tqcovus" {out}/1.blast.filter

# chrome locate
python3 {p_tools}/chrom_locate.py {out}/1.blast.filter -i {chrome_file} -o {out}/2.blast.chorme

# ps:去掉符合过滤要求，但是比对上染色体文件之外的blast记录
awk '$2!="None"' {out}/2.blast.chorme >{out}/3.blast.chorme.filter

# Rtab
python3 {p_tools}/filter_blast.py stat -f 1,2,15,3 {out}/3.blast.chorme.filter -o {out}/4.blast.same.Rtab

# show
chmod 744 {out}/4.blast.same.Rtab.stat

"""
    with open(f"{out}/run.sh", 'w', encoding='utf-8') as w:
        w.write(cmd)

    p = Popen(cmd, shell=True, encoding='utf-8', stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    print(stdout)
    open(f"{out}/run.sh.o", 'w', encoding='utf-8').write(stdout)
    open(f"{out}/run.sh.e", 'w', encoding='utf-8').write(stderr)


if __name__ == '__main__':
    main()
