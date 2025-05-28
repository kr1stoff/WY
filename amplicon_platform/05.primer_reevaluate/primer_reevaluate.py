#!/usr/bin/python3
# @auther: stone
# @date: 2023-07-18
import os
import click
import time
import logging
import glob
import subprocess


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


# == Option =====================================================================

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input', required=True, type=click.Path(), help="输入文件,table格式,分别为ID,")
@click.option('-t', '--p_type', default='nested', type=str, help="评估的类型, nested|tNGS|qPCR")
@click.option('-n', '--nt_database', default='/sdbb/bioinfor/yehui/nt/nt.fa', type=click.Path(), help="指定blast比对库")
@click.option('-d', '--database_path', required=False, type=click.Path(), help="/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge")
@click.option('-b', '--bacground_list', default='/sdbb/bioinfor/lanlei/script/amplicon_platform/src/all.bac.list.subtaxid', type=click.Path(), help="指定背景菌列表")
@click.option('--run/--no-run', default=False, show_default=True, help="是否直接运行")
@click.option('-o', '--outdir', default='./', type=click.Path(), help="输出目录,默认[./]")
def main(input, outdir, p_type, nt_database, database_path, bacground_list, run):
    input = os.path.abspath(input)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    Bin = os.path.dirname(os.path.abspath(__file__))

    shell_dir = f"{outdir}/00.shell"
    specificity_dir = f"{outdir}/01.specificity"
    inclusiveness_dir = f"{outdir}/02.inclusiveness"
    stat_dir = f"{outdir}/03.stat"

    mkdir_shell = f'mkdir -p {shell_dir} {specificity_dir} {inclusiveness_dir} {stat_dir} {outdir}/taxid_list'
    os.system(mkdir_shell)

    outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs staxid sscinames sstrand'"
    primer_blast_dir = os.path.abspath(f'{Bin}/../primer_blast')
    conservative_dir = os.path.abspath(f'{Bin}/../conservative')
    dimer_bin = os.path.abspath(f'{Bin}/../dimer')
    common_bin = os.path.abspath(f'{Bin}/../common_bin')
    database_path = '/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge'

    with open(input, 'r') as IN, open(f'{outdir}/outprimer.table.taxid', 'w') as T1, open(f'{outdir}/inprimer.table.taxid', 'w') as T2, open(f'{outdir}/probe.table.taxid', 'w') as T3, open(f'{outdir}/outprimer.fa', 'w') as F1, open(f'{outdir}/inprimer.fa', 'w') as F2, open(f'{outdir}/probe.fa', 'w') as F3:
        for line in IN:
            arr = line.strip().split('\t')
            name, in_f, in_r, probe, out_f, out_r, spe_taxid, inc_taxid = [
                arr[i] for i in [0, 1, 2, 3, 4, 5, 6, 7]]

            # 生成table
            print(*[f'{name}', f'{name}__F__out', f'{name}__R__out',
                  spe_taxid, inc_taxid], sep='\t', file=T1)
            print(*[f'{name}', f'{name}__F__in', f'{name}__R__in',
                  spe_taxid, inc_taxid], sep='\t', file=T2)
            print(*[f'{name}', f'{name}__P__in', f'{name}__N__N',
                  spe_taxid, inc_taxid], sep='\t', file=T3)

            # 生成fa
            print(f'>{name}__F__out\n{out_f}\n>{name}__R__out\n{out_r}', file=F1)
            print(f'>{name}__F__in\n{in_f}\n>{name}__R__in\n{in_r}', file=F2)
            print(f'>{name}__P__in\n{probe}', file=F3)

    with open(f'{shell_dir}/main.sh', 'w') as MAIN_shell:
        S1_shell_qPCR = f"""
######## inprimer #########
mkdir -p {specificity_dir}/inprimer
# cd {specificity_dir}/inprimer
python {Bin}/degenerate/IUPAC_translate.py -i {outdir}/inprimer.fa -o {outdir}/inprimer.degenerate.fa
python {common_bin}/blast_nt_cut.py --fasta {outdir}/inprimer.degenerate.fa --nt_database {nt_database} --out {specificity_dir}/inprimer/inprimer.blastn.out

perl {primer_blast_dir}/filter_blast.pl  {specificity_dir}/inprimer/inprimer.blastn.out > {specificity_dir}/inprimer/inprimer.blastn.out.filter
perl {primer_blast_dir}/split_workdir.pl {outdir}/inprimer.table.taxid {specificity_dir}/inprimer/inprimer.blastn.out.filter {specificity_dir}/inprimer

cut -f 4 {outdir}/inprimer.table.taxid | sort | uniq | while read i; do taxonkit list -i $i | awk '{{print $1}}' | grep -v '^$' > {outdir}/taxid_list/$i.taxid.list; done

cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/judge_primer_blast.pl {outdir}/taxid_list/$taxid_spe.taxid.list {specificity_dir}/inprimer/$id/$id.blast > {specificity_dir}/inprimer/$id/$id.blast.same; done
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/same_stat.pl {specificity_dir}/inprimer/$id/$id.blast.same > {specificity_dir}/inprimer/$id/$id.blast.same.stat; done
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/blast_primer_pair.pl {specificity_dir}/inprimer/$id/$id.input.table {specificity_dir}/inprimer/$id/$id.blast.same $id > {specificity_dir}/inprimer/$id/$id.out ; done
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/judge_primer_blast.pl {outdir}/taxid_list/$taxid_spe.taxid.list {specificity_dir}/inprimer/$id/$id.out 3 > {specificity_dir}/inprimer/$id/$id.out.same; done
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/same_stat.pl {specificity_dir}/inprimer/$id/$id.out.same > {specificity_dir}/inprimer/$id/$id.out.same.stat; done
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/percent_stat.pl $id {specificity_dir}/inprimer/$id/$id.input.table {specificity_dir}/inprimer/$id/$id.blast.same.stat {specificity_dir}/inprimer/$id/$id.out.same.stat > {specificity_dir}/inprimer/$id/$id.all.stat; done

######## probe #########
mkdir -p {specificity_dir}/probe
# cd {specificity_dir}/probe
python {Bin}/degenerate/IUPAC_translate.py -i {outdir}/probe.fa -o {outdir}/probe.degenerate.fa
python {common_bin}/blast_nt_cut.py --fasta {outdir}/probe.degenerate.fa --nt_database {nt_database} --out {specificity_dir}/probe/probe.blastn.out

perl {primer_blast_dir}/filter_blast.pl {specificity_dir}/probe/probe.blastn.out > {specificity_dir}/probe/probe.blastn.out.filter
perl {primer_blast_dir}/split_workdir.pl {outdir}/probe.table.taxid {specificity_dir}/probe/probe.blastn.out.filter {specificity_dir}/probe
cat {outdir}/probe.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/judge_primer_blast.pl {outdir}/taxid_list/$taxid_spe.taxid.list {specificity_dir}/probe/$id/$id.blast > {specificity_dir}/probe/$id/$id.blast.same; done
cat {outdir}/probe.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/same_stat.pl {specificity_dir}/probe/$id/$id.blast.same > {specificity_dir}/probe/$id/$id.blast.same.stat; done
cat {outdir}/probe.table.taxid | while read -r id f r taxid_spe taxid; do perl {conservative_dir}/probe_blast_stat.pl {specificity_dir}/probe/$id/$id.input.table {specificity_dir}/probe/$id/$id.blast.same.stat > {specificity_dir}/probe/$id/$id.all.stat ; done
"""
        S1_shell_tNGS = f"""
######## outprimer #########
mkdir -p {specificity_dir}/outprimer
# cd {specificity_dir}/outprimer
python {Bin}/degenerate/IUPAC_translate.py -i {outdir}/outprimer.fa -o {outdir}/outprimer.degenerate.fa
python {common_bin}/blast_nt_cut.py --fasta {outdir}/outprimer.degenerate.fa --nt_database {nt_database} --out {specificity_dir}/outprimer/outprimer.blastn.out

perl {primer_blast_dir}/filter_blast.pl  {specificity_dir}/outprimer/outprimer.blastn.out > {specificity_dir}/outprimer/outprimer.blastn.out.filter
perl {primer_blast_dir}/split_workdir.pl {outdir}/outprimer.table.taxid {specificity_dir}/outprimer/outprimer.blastn.out.filter {specificity_dir}/outprimer

cut -f 4 {outdir}/inprimer.table.taxid | sort | uniq | while read i; do taxonkit list -i $i | awk '{{print $1}}' | grep -v '^$' > {outdir}/taxid_list/$i.taxid.list; done

cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/judge_primer_blast.pl {outdir}/taxid_list/$taxid_spe.taxid.list {specificity_dir}/outprimer/$id/$id.blast > {specificity_dir}/outprimer/$id/$id.blast.same; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/same_stat.pl {specificity_dir}/outprimer/$id/$id.blast.same > {specificity_dir}/outprimer/$id/$id.blast.same.stat; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/blast_primer_pair.pl {specificity_dir}/outprimer/$id/$id.input.table {specificity_dir}/outprimer/$id/$id.blast.same $id > {specificity_dir}/outprimer/$id/$id.out ; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/judge_primer_blast.pl {outdir}/taxid_list/$taxid_spe.taxid.list {specificity_dir}/outprimer/$id/$id.out 3 > {specificity_dir}/outprimer/$id/$id.out.same; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/same_stat.pl {specificity_dir}/outprimer/$id/$id.out.same > {specificity_dir}/outprimer/$id/$id.out.same.stat; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/percent_stat.pl $id {specificity_dir}/outprimer/$id/$id.input.table {specificity_dir}/outprimer/$id/$id.blast.same.stat {specificity_dir}/outprimer/$id/$id.out.same.stat > {specificity_dir}/outprimer/$id/$id.all.stat; done
"""
        S1_shell_merge = f"""
######## merge #########
mkdir -p {specificity_dir}/merge
cut -f 1 {outdir}/outprimer.table.taxid | while read id; do python3 {Bin}/stat/merge_primer_blast.py --inprimer_blast {specificity_dir}/inprimer/$id/$id.out.same --outprimer_blast {specificity_dir}/outprimer/$id/$id.out.same --probe_blast {specificity_dir}/probe/$id/$id.blast.same > {specificity_dir}/merge/$id.merge.txt; done
cut -f 1 {outdir}/outprimer.table.taxid | while read id; do perl {primer_blast_dir}/same_stat.pl {specificity_dir}/merge/$id.merge.txt  > {specificity_dir}/merge/$id.specificity.txt; done
cut -f 1 {outdir}/outprimer.table.taxid | while read id; do python {Bin}/nonSpecific/deal_nonspecific.py {bacground_list} {specificity_dir}/merge/$id.merge.txt > {specificity_dir}/merge/$id.nonspecific.txt; done
"""
        S2_shell_qPCR = f"""
######## inprimer #########
mkdir -p {inclusiveness_dir}/inprimer
# cd {inclusiveness_dir}/inprimer

perl {conservative_dir}/split_workdir_blast_genome.pl {outdir}/inprimer.table.taxid {outdir}/inprimer.degenerate.fa {inclusiveness_dir}/inprimer
rm {inclusiveness_dir}/inprimer/inprimer.blastn.sh -rf
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do echo "blastn -task blastn-short -query {inclusiveness_dir}/inprimer/$id/$id.fa -db {database_path}/$taxid/all -num_threads 60 -word_size 9 -evalue 800000 -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt {outfmt} -out {inclusiveness_dir}/inprimer/$id/$id.blastn.out" >> {inclusiveness_dir}/inprimer/inprimer.blastn.sh; done

# perl {common_bin}/qsub_all_local.pl -b 1 -m 5 {inclusiveness_dir}/inprimer/inprimer.blastn.sh -d {inclusiveness_dir}/inprimer/qsub_{int(time.time())}
cat {inclusiveness_dir}/inprimer/inprimer.blastn.sh | parallel -j 5

cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/filter_blast.pl {inclusiveness_dir}/inprimer/$id/$id.blastn.out > {inclusiveness_dir}/inprimer/$id/$id.blastn.out.filter; done
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/blast_primer_pair.pl {inclusiveness_dir}/inprimer/$id/$id.input.table {inclusiveness_dir}/inprimer/$id/$id.blastn.out.filter $id > {inclusiveness_dir}/inprimer/$id/$id.out; done

cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {conservative_dir}/epcr_result_deal.pl {database_path}/$taxid/assembly_chromosome.tsv {inclusiveness_dir}/inprimer/$id/$id.out {inclusiveness_dir}/inprimer/$id/$id.assembly_chromosome.copy {inclusiveness_dir}/inprimer/$id/$id.out.chr {inclusiveness_dir}/inprimer/$id/$id.all.stat; done

######## probe #########
mkdir -p {inclusiveness_dir}/probe
# cd {inclusiveness_dir}/probe
rm {inclusiveness_dir}/probe/probe.blastn.sh -rf
perl {conservative_dir}/split_workdir_blast_genome_probe.pl {outdir}/probe.table.taxid {outdir}/probe.degenerate.fa {inclusiveness_dir}/probe
cat {outdir}/inprimer.table.taxid | while read -r id f r taxid_spe taxid; do echo "blastn -task blastn-short -query {inclusiveness_dir}/probe/$id/$id.fa -db {database_path}/$taxid/all -num_threads 60 -word_size 9 -evalue 800000 -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt {outfmt} -out {inclusiveness_dir}/probe/$id/$id.blastn.out" >> {inclusiveness_dir}/probe/probe.blastn.sh; done

# perl {common_bin}/qsub_all_local.pl -b 1 -m 5 {inclusiveness_dir}/probe/probe.blastn.sh -d {inclusiveness_dir}/probe/qsub_{int(time.time())}
cat {inclusiveness_dir}/probe/probe.blastn.sh | parallel -j 5

cat {outdir}/probe.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/filter_blast.pl {inclusiveness_dir}/probe/$id/$id.blastn.out > {inclusiveness_dir}/probe/$id/$id.blastn.out.filter; done

cat {outdir}/probe.table.taxid | while read -r id f r taxid_spe taxid; do perl {conservative_dir}/epcr_result_deal.pl {database_path}/$taxid/assembly_chromosome.tsv {inclusiveness_dir}/probe/$id/$id.blastn.out.filter {inclusiveness_dir}/probe/$id/$id.assembly_chromosome.copy {inclusiveness_dir}/probe/$id/$id.blast.chr {inclusiveness_dir}/probe/$id/$id.blast.stat; done

cat {outdir}/probe.table.taxid | while read -r id f r taxid_spe taxid; do perl {conservative_dir}/probe_blast_stat.pl {inclusiveness_dir}/probe/$id/$id.input.table {inclusiveness_dir}/probe/$id/$id.blast.stat > {inclusiveness_dir}/probe/$id/$id.all.stat ; done
"""
        S2_shell_tNGS = f"""
######## outprimer #########
mkdir -p {inclusiveness_dir}/outprimer
# cd {inclusiveness_dir}/outprimer
perl {conservative_dir}/split_workdir_blast_genome.pl {outdir}/outprimer.table.taxid {outdir}/outprimer.degenerate.fa {inclusiveness_dir}/outprimer
rm {inclusiveness_dir}/outprimer/outprimer.blastn.sh -rf
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do echo "blastn -task blastn-short -query {inclusiveness_dir}/outprimer/$id/$id.fa -db {database_path}/$taxid/all -num_threads 80 -word_size 9 -evalue 800000 -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt {outfmt} -out {inclusiveness_dir}/outprimer/$id/$id.blastn.out" >> {inclusiveness_dir}/outprimer/outprimer.blastn.sh; done

# perl {common_bin}/qsub_all_local.pl -b 1 -m 5 {inclusiveness_dir}/outprimer/outprimer.blastn.sh -d {inclusiveness_dir}/outprimer/qsub_{int(time.time())}
cat {inclusiveness_dir}/outprimer/outprimer.blastn.sh | parallel -j 5

cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/filter_blast.pl {inclusiveness_dir}/outprimer/$id/$id.blastn.out > {inclusiveness_dir}/outprimer/$id/$id.blastn.out.filter; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {primer_blast_dir}/blast_primer_pair.pl {inclusiveness_dir}/outprimer/$id/$id.input.table {inclusiveness_dir}/outprimer/$id/$id.blastn.out.filter $id > {inclusiveness_dir}/outprimer/$id/$id.out; done

cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do perl {conservative_dir}/epcr_result_deal.pl {database_path}/$taxid/assembly_chromosome.tsv {inclusiveness_dir}/outprimer/$id/$id.out {inclusiveness_dir}/outprimer/$id/$id.assembly_chromosome.copy {inclusiveness_dir}/outprimer/$id/$id.out.chr {inclusiveness_dir}/outprimer/$id/$id.all.stat; done
"""
        S2_shell_merge = f"""
######## merge #########
mkdir -p {inclusiveness_dir}/merge
cut -f 1 {outdir}/outprimer.table.taxid | while read id;do python3 {Bin}/stat/epcr_inclusiveness_stat.py -i $id -1 {inclusiveness_dir}/inprimer/$id/$id.out.chr -2 {inclusiveness_dir}/outprimer/$id/$id.out.chr -3 {inclusiveness_dir}/probe/$id/$id.blast.chr -c {inclusiveness_dir}/probe/$id/$id.assembly_chromosome.copy > {inclusiveness_dir}/merge/$id.inclusiveness.txt; done
"""

        S3_shell_nested = f"""
# cd stat_dir
rm {stat_dir}/*all.stat {stat_dir}/*stat.txt -rf
# specificity
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/01.specificity/inprimer -o {stat_dir}/inprimer.specificity.all.stat
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/01.specificity/probe -o {stat_dir}/probe.specificity.all.stat
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/01.specificity/outprimer -o {stat_dir}/outprimer.specificity.all.stat
# inclusiveness
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/02.inclusiveness/inprimer -o {stat_dir}/inprimer.inclusiveness.all.stat
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/02.inclusiveness/probe -o {stat_dir}/probe.inclusiveness.all.stat
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/02.inclusiveness/outprimer -o {stat_dir}/outprimer.inclusiveness.all.stat
# merge
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/01.specificity/merge -o {stat_dir}/merge.specificity.stat.txt -t '$id.specificity.txt'
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/02.inclusiveness/merge -o {stat_dir}/merge.inclusiveness.stat.txt -t '$id.inclusiveness.txt'
# merge_final
python {Bin}/stat/merge_evaluate_result.py nested -i {input} -o {stat_dir}/primer.evaluate.txt -1 {stat_dir}/inprimer.specificity.all.stat -2 {stat_dir}/probe.specificity.all.stat -3 {stat_dir}/outprimer.specificity.all.stat -4 {stat_dir}/merge.specificity.stat.txt -5 {stat_dir}/inprimer.inclusiveness.all.stat -6 {stat_dir}/probe.inclusiveness.all.stat -7 {stat_dir}/outprimer.inclusiveness.all.stat -8 {stat_dir}/merge.inclusiveness.stat.txt
# nonspecific
python {Bin}/stat/cat_all_stat.py -i {outdir}/outprimer.table.taxid -f {outdir}/01.specificity/merge -o {stat_dir}/merge.nonspecific.all.stat -t '$id.nonspecific.txt'
# dimer
mkdir -p {stat_dir}/dimer
python {Bin}/degenerate/IUPAC_translate.py --add -i {outdir}/inprimer.fa -o {stat_dir}/dimer/inprimer.fa
python {Bin}/degenerate/IUPAC_translate.py --add -i {outdir}/outprimer.fa -o {stat_dir}/dimer/outprimer.fa
python {dimer_bin}/DimerAnalysis.py -i {stat_dir}/dimer/inprimer.fa -o {stat_dir}/dimer/
python {dimer_bin}/DimerAnalysis.py -i {stat_dir}/dimer/outprimer.fa -o {stat_dir}/dimer/
python {dimer_bin}/IUPAC_result_merge.py -i {stat_dir}/dimer/inprimer.mfeprimer.filter -o {stat_dir}/dimer/inprimer.mfeprimer.filter.merge
python {dimer_bin}/IUPAC_result_merge.py -i {stat_dir}/dimer/outprimer.mfeprimer.filter -o {stat_dir}/dimer/outprimer.mfeprimer.filter.merge
python {dimer_bin}/dimer_result_stat.py -i {stat_dir}/dimer/inprimer.mfeprimer.filter.merge -t {outdir}/inprimer.table.taxid -s in -o {stat_dir}/dimer/inprimer.dimer.txt
python {dimer_bin}/dimer_result_stat.py -i {stat_dir}/dimer/outprimer.mfeprimer.filter.merge -t {outdir}/outprimer.table.taxid -s out -o {stat_dir}/dimer/outprimer.dimer.txt
python {Bin}/stat/add_nonspecific_dimer.py --stat {stat_dir}/primer.evaluate.txt --nonspecific {stat_dir}/merge.nonspecific.all.stat --dimer_in {stat_dir}/dimer/inprimer.dimer.txt --dimer_out {stat_dir}/dimer/outprimer.dimer.txt --out {stat_dir}/primer.merge.txt
"""

        S3_shell_qPCR = f"""
# cd {stat_dir}
rm {stat_dir}/*all.stat {stat_dir}/*stat.txt -rf
# specificity
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do cat {outdir}/01.specificity/inprimer/$id/$id.all.stat | grep -v '#' >> {stat_dir}/inprimer.specificity.all.stat; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do cat {outdir}/01.specificity/probe/$id/$id.all.stat | grep -v '#' >> {stat_dir}/probe.specificity.all.stat; done
# inclusiveness
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do cat {outdir}/02.inclusiveness/inprimer/$id/$id.all.stat | grep -v '#' >> {stat_dir}/inprimer.inclusiveness.all.stat; done
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do cat {outdir}/02.inclusiveness/probe/$id/$id.all.stat | grep -v '#' >> {stat_dir}/probe.inclusiveness.all.stat; done
# merge_final
python {Bin}/stat/merge_evaluate_result.py qpcr -i {input} -o {stat_dir}/primer.merge.txt -1 {stat_dir}/inprimer.specificity.all.stat -2 {stat_dir}/probe.specificity.all.stat -3 {stat_dir}/inprimer.inclusiveness.all.stat -4 {stat_dir}/probe.inclusiveness.all.stat
"""

        S3_shell_tNGS = f"""
# cd {stat_dir}
rm {stat_dir}/*all.stat {stat_dir}/*stat.txt -rf
# specificity
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do cat {outdir}/01.specificity/outprimer/$id/$id.all.stat | grep -v '#' >> {stat_dir}/outprimer.specificity.all.stat; done
# inclusiveness
cat {outdir}/outprimer.table.taxid | while read -r id f r taxid_spe taxid; do cat {outdir}/02.inclusiveness/outprimer/$id/$id.all.stat | grep -v '#' >> {stat_dir}/outprimer.inclusiveness.all.stat; done
# merge_final
python {Bin}/stat/merge_evaluate_result.py tngs -i {input} -o {stat_dir}/primer.merge.txt -1 {stat_dir}/outprimer.specificity.all.stat -2 {stat_dir}/outprimer.inclusiveness.all.stat
"""

        shells = []
        if (p_type == "nested"):
            generate_shell(f'{shell_dir}/S1.specificity_qPCR.sh',
                           S1_shell_qPCR, shells)
            generate_shell(f'{shell_dir}/S1.specificity_tNGS.sh',
                           S1_shell_tNGS, shells)
            generate_shell(f'{shell_dir}/S1.specificity_merge.sh',
                           S1_shell_merge, shells)
            generate_shell(f'{shell_dir}/S2.inclusiveness_qPCR.sh',
                           S2_shell_qPCR, shells)
            generate_shell(f'{shell_dir}/S2.inclusiveness_tNGS.sh',
                           S2_shell_tNGS, shells)
            generate_shell(f'{shell_dir}/S2.inclusiveness_merge.sh',
                           S2_shell_merge, shells)
            generate_shell(f'{shell_dir}/S3.stat.sh', S3_shell_nested, shells)
        elif (p_type == "qPCR"):
            generate_shell(f'{shell_dir}/S1.specificity_qPCR.sh',
                           S1_shell_qPCR, shells)
            generate_shell(f'{shell_dir}/S2.inclusiveness_qPCR.sh',
                           S2_shell_qPCR, shells)
            generate_shell(f'{shell_dir}/S3.stat.sh', S3_shell_qPCR, shells)
        elif (p_type == "tNGS"):
            generate_shell(f'{shell_dir}/S1.specificity_tNGS.sh',
                           S1_shell_tNGS, shells)
            generate_shell(f'{shell_dir}/S2.inclusiveness_tNGS.sh',
                           S2_shell_tNGS, shells)
            generate_shell(f'{shell_dir}/S3.stat.sh', S3_shell_tNGS, shells)
        else:
            logging.error(f"Invalid type: {type}")
        main_shell = " && \\\n".join(shells)
        print(main_shell, file=MAIN_shell)
    if run:
        try:
            subprocess.run(
                f"bash {shell_dir}/main.sh 1>{shell_dir}/main.sh.e 2>{shell_dir}/main.sh.o", shell=True)
        except subprocess.CalledProcessError as e:
            print("命令执行失败:", e)


if __name__ == '__main__':
    main()
