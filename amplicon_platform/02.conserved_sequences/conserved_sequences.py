# -*- coding: utf-8 -*-
# @Author: Stone
# @Date:   2023-08-08 16:14

import logging
import os
import yaml
import click
import glob
import subprocess

#### Some Functions ####
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
__version__ = '1.0.0'

########################

# === Common Fuction =============================================================================

Bin = os.path.dirname(os.path.abspath(__file__))
common_bin = os.path.abspath(f'{Bin}/../common_bin')


def generate_shell(output_shell, content, outshell, finish_string=None):
    if finish_string is None:
        finish_string = "Still_waters_run_deep"

    # Remove existing output_shell files
    existing_files = glob.glob(f"{output_shell}.*")
    for file in existing_files:
        os.remove(file)

    # Write content to the output_shell file
    with open(output_shell, "w") as out_file:
        out_file.write("#!/bin/bash\n")
        out_file.write("echo ==========start at : `date` ==========\n")
        out_file.write("set -e \n")
        out_file.write(f"{content} \n")
        out_file.write("echo ==========end at : `date` ========== && \n")
        out_file.write(f"echo {finish_string} 1>&2 && \\\n")
        out_file.write(f"echo {finish_string} > {output_shell}.sign\n")

    run_command = f"bash {output_shell} 1>{output_shell}.e 2>{output_shell}.o"
    outshell.append(run_command)


def inclusiveness(fa, reference_dir, outdir, inc_num):
    head_num = inc_num + 1
    outfmt_inc = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus'"

    inclusiveness_dir = f"{outdir}/02.inclusiveness"

    S2_shell = f"""
# blast
blastn -db {reference_dir}/all -evalue 1e-5 -num_threads 40 -max_target_seqs 5000 -outfmt {outfmt_inc} -query {fa} -out {inclusiveness_dir}/0.blast

# filter --cov 95 --identity 95
python3 {Bin}/inclusiveness/filter_blast.py filter  {inclusiveness_dir}/0.blast -f 15,3 --cov 95 --identity 95 -o {inclusiveness_dir}/1.blast.filter

# add-header
# sed -i "1iqaccver\\tsaccver\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\tqcovs\\tqcovhsp\\tqcovus" {inclusiveness_dir}/1.blast.filter

# chrome locate
python3 {Bin}/inclusiveness/chrom_locate.py {inclusiveness_dir}/1.blast.filter -i {reference_dir}/assembly_chromosome.tsv -o {inclusiveness_dir}/2.blast.chorme

# ps:去掉符合过滤要求,但是比对上染色体文件之外的blast记录
awk '$2!="None"' {inclusiveness_dir}/2.blast.chorme > {inclusiveness_dir}/3.blast.chorme.filter

# Rtab
python3 {Bin}/inclusiveness/filter_blast.py stat -f 1,2,15,3 {inclusiveness_dir}/3.blast.chorme.filter -o {inclusiveness_dir}/4.blast.same.Rtab

# filter
awk '$4>0.8' {inclusiveness_dir}/4.blast.same.Rtab.stat | head -n {head_num} > {inclusiveness_dir}/4.blast.same.Rtab.stat.filter

# get filter fa
perl {common_bin}/fishInWinter.pl -ff fasta {inclusiveness_dir}/4.blast.same.Rtab.stat.filter {fa} > {inclusiveness_dir}/inclusiveness.filter.fa
"""
    return S2_shell


def specificity(fa, outdir, taxid, nt_database):

    specificity_dir = f"{outdir}/03.specificity"

    S3_shell = f"""
# 生成某个taxid下面所有的taxid,即taxid.list文件,后面要用
taxonkit list -i {taxid} | awk '{{print $1}}' | grep -v '^$' > {specificity_dir}/taxid.list
# 比对到nt库
python {common_bin}/blast_nt_cut.py --fasta {fa} --nt_database {nt_database} --long --out {specificity_dir}/blastn.out
# step1过滤,80%的coverage和90%的identity
perl {Bin}/specificity/step1.filter.cov80_and_iden90.pl {specificity_dir}/blastn.out > {specificity_dir}/blastn.out.step1
# step2,加上same和diff
perl {Bin}/specificity/step2.judge_spe.pl {specificity_dir}/taxid.list {specificity_dir}/blastn.out.step1 > {specificity_dir}/blastn.out.step2
# step3,统计same和diff
perl {Bin}/specificity/step3.filter_consensus.pl {specificity_dir}/blastn.out.step2 > {specificity_dir}/blastn.out.step3
# percentage>80%
cat {specificity_dir}/blastn.out.step3 | awk '$2>80' > {specificity_dir}/blastn.out.step3.filter
"""
    return S3_shell


def merge(spe_result, inc_result, outdir, taxid, start_num, name2taxid):

    merge_dir = f"{outdir}/04.merge"
    S4_shell = f"""
# 合并特异性，核心基因结果，第二遍包容性评估结果
perl {Bin}/merge/merge_specific_conserved.pl {spe_result} {inc_result} > {merge_dir}/result.txt

# 根据结果数手动调整参数阈值
python {Bin}/merge/resion_level.py --input_result {merge_dir}/result.txt --output_result {merge_dir}/result.txt.level
python {Bin}/merge/get_filter_result_bed.py --input_result {merge_dir}/result.txt.level --output_result {merge_dir}/result.txt.level.filter --start_num {start_num} --config {name2taxid} --taxid {taxid}
"""
    return S4_shell


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    """
    Pipeline for conserver sequence idetity
    """
    pass

# ===== *gene ===========================================================================================================


@click.command()
@click.option('-r', '--reference_dir',
              type=click.Path(),
              required=True,
              help='The ref fasta file, not fasta')
@click.option('-g', '--gff',
              type=click.Path(),
              required=True,
              help='The gene gff file')
@click.option('-t', '--taxid',
              type=str,
              required=True,
              help='The species taxid')
@click.option('-o', '--outdir',
              type=click.Path(),
              default="./",
              help='The outdir')
@click.option('-w', '--windows',
              type=int,
              default=300,
              help='The sliding window method: windows size')
@click.option('-s', '--step',
              type=int,
              default=300,
              help='The sliding window method: step size')
@click.option('-n', '--inc_num',
              type=int,
              default=10000,
              help='The resion number for specificity analysis')
@click.option('-st', '--start_num',
              type=int,
              default=1,
              help='The bed name start id')
@click.option('-d', '--nt_database',
              type=click.Path(),
              help='NT database path, default: yaml file')
@click.option("--run/--no-run",
              default=False,
              show_default=True,
              help="Whether run the script directly")
@click.option('-y', '--f_yaml',
              type=click.Path(),
              required=True,
              help='The config yaml file')
def gene(reference_dir, gff, taxid, outdir, windows, step, inc_num, start_num, nt_database, run, f_yaml):
    """
    Choose gene to get conserved sequnence
    """
    info = yaml.safe_load(open(f_yaml, 'r'))
    if nt_database:
        info['nt_database'] = nt_database
    logging.info("Choose gene to get conserved sequnence")
    reference_dir = os.path.abspath(reference_dir)
    gff = os.path.abspath(gff)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    shell_dir = f"{outdir}/00.shell"
    sliding_dir = f"{outdir}/01.sliding"
    inclusiveness_dir = f"{outdir}/02.inclusiveness"
    specificity_dir = f"{outdir}/03.specificity"
    merge_dir = f"{outdir}/04.merge"

    mkdir_shell = f'mkdir -p {shell_dir} {sliding_dir} {specificity_dir} {inclusiveness_dir} {merge_dir}'
    os.system(mkdir_shell)

    S1_shell = f"""
# gene的bed文件
awk '$3=="CDS"' {gff} | awk -F '\\t' -v 'OFS=\\t' '$5-$4>200{{print $1,$4-1,$5}}' > {sliding_dir}/ref_gene.bed
# 生成基因区的fa
bedtools getfasta -fi {reference_dir}/ref.fna -fo {sliding_dir}/ref_gene.fa -bed {sliding_dir}/ref_gene.bed
# 滑窗切fa
seqkit sliding -j 8 -s {step} -W {windows} {sliding_dir}/ref_gene.fa -o {sliding_dir}/ref_gene.sliding.fna
"""
    S2_shell = inclusiveness(
        f'{sliding_dir}/ref_gene.sliding.fna', reference_dir, outdir, inc_num)
    S3_shell = specificity(
        f'{inclusiveness_dir}/inclusiveness.filter.fa', outdir, taxid, info['nt_database'])
    S4_shell = merge(
        f'{specificity_dir}/blastn.out.step3.filter',
        f'{inclusiveness_dir}/4.blast.same.Rtab.stat.filter', outdir, taxid, start_num, info['name2taxid'])

    shells = []
    generate_shell(f'{shell_dir}/S1.sliding.sh',
                   S1_shell, shells)
    generate_shell(f'{shell_dir}/S2.inclusiveness.sh',
                   S2_shell, shells)
    generate_shell(f'{shell_dir}/S3.specificity.sh',
                   S3_shell, shells)
    generate_shell(f'{shell_dir}/S4.merge.sh',
                   S4_shell, shells)
    main_shell = " && \\\n".join(shells)
    with open(f'{shell_dir}/main.sh', 'w') as MAIN_shell:
        print(main_shell, file=MAIN_shell)
    if run:
        subprocess.run(
            f"bash {shell_dir}/main.sh 1>{shell_dir}/main.sh.e 2>{shell_dir}/main.sh.o", shell=True)


# ===== *genome ===========================================================================================================
@click.command()
@click.option('-d', '--reference_dir',
              type=click.Path(),
              required=True,
              help='The ref fasta file, not fasta')
@click.option('-t', '--taxid',
              type=str,
              required=True,
              help='The species taxid')
@click.option('-o', '--outdir',
              type=click.Path(),
              default="./",
              help='The outdir')
@click.option('-w', '--windows',
              type=int,
              default=300,
              help='The sliding window method: windows size')
@click.option('-s', '--step',
              type=int,
              default=300,
              help='The sliding window method: step size')
@click.option('-n', '--inc_num',
              type=int,
              default=10000,
              help='The resion number for specificity analysis')
@click.option('-st', '--start_num',
              type=int,
              default=1,
              help='The bed name start id, default: 1')
@click.option('-d', '--nt_database',
              type=click.Path(),
              help='NT database path, default: yaml file')
@click.option("--run/--no-run",
              default=False,
              show_default=True,
              help="Whether run the script directly")
@click.option('-y', '--f_yaml',
              type=click.Path(),
              default=f'{Bin}/conserved_sequences.yaml',
              help='The config yaml file')
def genome(reference_dir, taxid, outdir, windows, step, inc_num, nt_database, start_num, run, f_yaml):
    """
    Choose gene to get conserved sequnence
    """
    info = yaml.safe_load(open(f_yaml, 'r'))
    if nt_database:
        info['nt_database'] = nt_database
    logging.info("Choose genome to get conserved sequnence")
    reference_dir = os.path.abspath(reference_dir)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    shell_dir = f"{outdir}/00.shell"
    sliding_dir = f"{outdir}/01.sliding"
    inclusiveness_dir = f"{outdir}/02.inclusiveness"
    specificity_dir = f"{outdir}/03.specificity"
    merge_dir = f"{outdir}/04.merge"

    mkdir_shell = f'mkdir -p {shell_dir} {sliding_dir} {specificity_dir} {inclusiveness_dir} {merge_dir}'
    os.system(mkdir_shell)

    S1_shell = f"""
# 滑窗切fa
seqkit sliding -j 8 -s {step} -W {windows} {reference_dir}/ref.fna -o {sliding_dir}/ref_genome.sliding.fna
"""
    S2_shell = inclusiveness(
        f'{sliding_dir}/ref_genome.sliding.fna', reference_dir, outdir, inc_num)
    S3_shell = specificity(
        f'{inclusiveness_dir}/inclusiveness.filter.fa', outdir, taxid, info['nt_database'])
    S4_shell = merge(
        f'{specificity_dir}/blastn.out.step3.filter',
        f'{inclusiveness_dir}/4.blast.same.Rtab.stat.filter', outdir, taxid, start_num, info['name2taxid'])

    shells = []
    generate_shell(f'{shell_dir}/S1.sliding.sh',
                   S1_shell, shells)
    generate_shell(f'{shell_dir}/S2.inclusiveness.sh',
                   S2_shell, shells)
    generate_shell(f'{shell_dir}/S3.specificity.sh',
                   S3_shell, shells)
    generate_shell(f'{shell_dir}/S4.merge.sh',
                   S4_shell, shells)
    main_shell = " && \\\n".join(shells)
    with open(f'{shell_dir}/main.sh', 'w') as MAIN_shell:
        print(main_shell, file=MAIN_shell)
    if run:
        subprocess.run(
            f"bash {shell_dir}/main.sh 1>{shell_dir}/main.sh.e 2>{shell_dir}/main.sh.o", shell=True)


cli.add_command(gene)
cli.add_command(genome)

if __name__ == "__main__":
    cli()
