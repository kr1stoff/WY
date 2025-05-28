# -*- coding: utf-8 -*-
# @Author: Stone
# @Date:   2023-08-08 16:14

import logging
import os
import yaml
import click

#### Some Functions ####
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
__version__ = '1.0.0'

########################

# === Common Fuction =============================================================================

Bin = os.path.dirname(os.path.abspath(__file__))

###################################################################################################
# 数据库下载
###################################################################################################


def Download(taxid, ptype, outdir):
    # logging.info("Start download genome")
    download_shell = f"""
python {Bin}/01.download/bin/downloader.py --taxonid {taxid} --ptype {ptype} --db {outdir} --name {taxid} --run
python {Bin}/01.download/bin/builder.py --taxonid {taxid} --db {outdir} --name {taxid}
"""
    return download_shell

###################################################################################################
# 保守区域预测
###################################################################################################


def Conserved_sequences(ptype, reference_dir, gff, outdir, taxid, windows, step, inc_num, nt_database, f_yaml):
    # logging.info("Start finding conserved sequnence")
    if (ptype == 'genome'):
        Conserved_sequences_shell = f"""
python {Bin}/02.conserved_sequences/conserved_sequences.py genome --reference_dir {reference_dir} --taxid {taxid} --outdir {outdir} --windows {windows} --step {step} --inc_num {inc_num} --nt_database {nt_database} --f_yaml {f_yaml} --run
"""
    elif (ptype == 'gene'):
        Conserved_sequences_shell = f"""
python {Bin}/02.conserved_sequences/conserved_sequences.py gene --reference_dir {reference_dir} --gff {gff} --taxid {taxid} --outdir {outdir} --windows {windows} --step {step} --inc_num {inc_num} --nt_database {nt_database} --f_yaml {f_yaml} --run
"""
    else:
        logging.warning("没有这种类型")

    return Conserved_sequences_shell

###################################################################################################
# 引物设计和评估大模块
###################################################################################################
# 内引物设计


def qPcr_Design(reference, bed, outdir):
    # logging.info("Start qPcr primer design")
    qPcr_Design_shell = f"""
python {Bin}/03.primer_evaluate/qpcrDesign/main.py qpcr-batch --reference {reference} --bed {bed} --out {outdir}
"""
    return qPcr_Design_shell

# 内引物评估


def qPcr_Evaluate(primer_file, genome_dir, name2taxid, nt_database, evalue, outdir):
    qPcr_Evaluate_shell = f"""
perl {Bin}/03.primer_evaluate/primer_evaluate_qpcr.pl -primer {primer_file} -genome_dir {genome_dir} -name2taxid {name2taxid} -nt_database {nt_database} -evalue {evalue} -outdir {outdir}
"""
    return qPcr_Evaluate_shell

# 外引物设计


def tNGS_Design(tNGS_bed, outdir):
    tNGS_Design_shell = f"""
perl {Bin}/03.primer_evaluate/primer_design.pl -regions {tNGS_bed}  -outdir {outdir}
"""
    return tNGS_Design_shell

# 外引物评估


def tNGS_Evaluate(primer_file, reference_dir, name2taxid, nt_database, evalue, outdir):
    tNGS_Evaluate_shell = f"""
perl {Bin}/03.primer_evaluate/primer_evaluate_pcr.pl -primer {primer_file} -genome_dir {reference_dir} -name2taxid {name2taxid} -nt_database {nt_database} -evalue {evalue} -outdir {outdir}
"""
    return tNGS_Evaluate_shell

# 内引物和外引物结果合并


def Merge_qPCR_tNGS(chinese, qPCR_xls, tNGS_xls, ref, base_num, outdir):
    Merge_qPCR_tNGS_shell = f"""
perl {Bin}/03.primer_evaluate/stat/merge_primer_result.pl -n {chinese} -i {qPCR_xls} -o {tNGS_xls} -r {ref} -l {base_num} -d {outdir}
"""
    return Merge_qPCR_tNGS_shell

###################################################################################################
# 巢式PCR模块
# 后续使用经常是内引物设计+评估，外引物设计+评估，都一起跑，所以合并成巢式PCR的脚本
###################################################################################################


def Nested_PCR(reference_dir, bed, outdir, name2taxid, nt_database, evalue, chinese, base_num):
    qPcr_Design_shell = qPcr_Design(
        f'{reference_dir}/ref.fna', bed, f'{outdir}/qPcr_Design')
    qPcr_Evaluate_shell = qPcr_Evaluate(
        f'{outdir}/qPcr_Design/exp.txt', reference_dir, name2taxid, nt_database, evalue, f'{outdir}/qPcr_Evaluate')
    tNGS_Design_shell = tNGS_Design(
        f'{outdir}/qPcr_Evaluate/03.stat/2tNGS.bed', f'{outdir}/tNGS_Design')
    tNGS_Evaluate_shell = tNGS_Evaluate(
        f'{outdir}/tNGS_Design/2tNGS_primers_combination_1_info.txt', reference_dir, name2taxid, nt_database, evalue, f'{outdir}/tNGS_Evaluate')
    Merge_qPCR_tNGS_shell = Merge_qPCR_tNGS(
        chinese, f'{outdir}/qPcr_Evaluate/03.stat/qPCR.primer.check.xls', f'{outdir}/tNGS_Evaluate/03.stat/out.primer.check.xls', f'{reference_dir}/ref.fna', base_num, f'{outdir}/Merge')
    Nested_PCR_shell = qPcr_Design_shell + qPcr_Evaluate_shell + \
        tNGS_Design_shell + tNGS_Evaluate_shell + Merge_qPCR_tNGS_shell
    return Nested_PCR_shell

###################################################################################################
# 简并碱基处理
###################################################################################################


def Degenerate_Bases(merge, allfna, outdir):
    merge_name = os.path.basename(merge)
    merge_name_prefix = os.path.splitext(merge_name)[0]
    Degenerate_Bases_shell = f"""
python {Bin}/04.degenerate_base/main.py translate-degenerate-base --primertable {merge} --allfna {allfna} --outdir {outdir}
python {Bin}/04.degenerate_base/xls2reevaluate_table.py {outdir}/{merge_name_prefix}.degenerate.xls > {outdir}/primer2taxid.table
"""
    return Degenerate_Bases_shell

###################################################################################################
# 引物重新评估
###################################################################################################


def Primer_Reevaluate(input_table, outdir, nt_database, database_path):
    Primer_Reevaluate_shell = f"""
python {Bin}/05.primer_reevaluate/primer_reevaluate.py --input {input_table} --outdir {outdir} --nt_database {nt_database} --database_path {database_path} --run
"""
    return Primer_Reevaluate_shell


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    """
    Pipeline for amplicon platform
    """
    pass

# ===== *gene ===========================================================================================================


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-y', '--f_yaml',
              type=click.Path(),
              required=True,
              help='The config yaml file')
@click.option('-o', '--outdir',
              type=click.Path(),
              default="./",
              help='The outdir')
def Amplicon(f_yaml, outdir):
    """
    Choose gene to get conserved sequnence
    """
    step2cmd = {'1': 'Download',
                '2': 'Conserved_sequences',
                '3': 'Nested_PCR',
                '4': 'Degenerate_Bases',
                '5': 'Primer_Reevaluate',
                '6': 'Primer_select'}
    f_yaml = os.path.abspath(f_yaml)
    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    # 读yaml文件
    info = yaml.safe_load(open(f_yaml, 'r'))
    os.system(f'mkdir -p {outdir}/00.shell')
    # 读run，看看跑哪几步
    run_yaml = info['run']
    print(f'Step: {run_yaml}')
    steps = list(str(run_yaml))
    # print(steps)
    first_step = steps[0]
    last_step = steps[-1]
    for step in steps:
        cmd_shell = step2cmd[step]
        mkdir_shell = f'mkdir -p 0{step}.{cmd_shell}'
        os.system(mkdir_shell)
    if (first_step == '1'):
        outshells = []
        taxid = info['taxid']
        # reference_dir = os.path.abspath(info['reference_dir'])
        reference_dir = f'{outdir}/01.Download/{taxid}'

        Download_shell = Download(
            info['taxid'], info['ptype'], f'{outdir}/01.Download')
        outshells.append(Download_shell)

        if (int(last_step) >= 2):
            Conserved_sequences_shell = Conserved_sequences(
                info['resion_type'], reference_dir, info['gff'], f'{outdir}/02.Conserved_sequences', info['taxid'], info['windows'], info['step'], info['inc_num'], info['nt_database'], f_yaml)
            outshells.append(Conserved_sequences_shell)

        if (int(last_step) >= 3):
            Nested_PCR_shell = Nested_PCR(reference_dir, f'{outdir}/02.Conserved_sequences/04.merge/result.txt.level.filter',
                                          f'{outdir}/03.Nested_PCR', info['name2taxid'], info['nt_database'], info['evalue'],
                                          info['chinese'], info['base_num'])
            outshells.append(Nested_PCR_shell)

        if (int(last_step) >= 4):
            chinese = info['chinese']
            Degenerate_Bases_shell = Degenerate_Bases(
                f'{outdir}/03.Nested_PCR/Merge/{chinese}_merge_primer.xls', f'{reference_dir}/all.fna', f'{outdir}/04.Degenerate_Bases/')
            outshells.append(Degenerate_Bases_shell)

        if (int(last_step) >= 5):
            Primer_Reevaluate_shell = Primer_Reevaluate(
                f'{outdir}/04.Degenerate_Bases/primer2taxid.table', f'{outdir}/05.Primer_Reevaluate', info['nt_database'], info['database_path'])
            outshells.append(Primer_Reevaluate_shell)

        with open(f'{outdir}/00.shell/run{run_yaml}.sh', 'w') as SHELL1:
            run_shell = "\n".join(outshells)
            print(run_shell, file=SHELL1)

    if (first_step == '2'):
        outshells = []
        reference_dir = os.path.abspath(info['reference_dir'])

        Conserved_sequences_shell = Conserved_sequences(
            info['resion_type'], reference_dir, info['gff'], f'{outdir}/02.Conserved_sequences', info['taxid'], info['windows'], info['step'], info['inc_num'], info['nt_database'], f_yaml)
        outshells.append(Conserved_sequences_shell)

        if (int(last_step) >= 3):
            Nested_PCR_shell = Nested_PCR(reference_dir, f'{outdir}/02.Conserved_sequences/04.merge/result.txt.level.filter',
                                          f'{outdir}/03.Nested_PCR', info['name2taxid'], info['nt_database'], info['evalue'],
                                          info['chinese'], info['base_num'])
            outshells.append(Nested_PCR_shell)

        if (int(last_step) >= 4):
            chinese = info['chinese']
            Degenerate_Bases_shell = Degenerate_Bases(
                f'{outdir}/03.Nested_PCR/Merge/{chinese}_merge_primer.xls', f'{reference_dir}/all.fna', f'{outdir}/04.Degenerate_Bases/')
            outshells.append(Degenerate_Bases_shell)

        if (int(last_step) >= 5):
            Primer_Reevaluate_shell = Primer_Reevaluate(
                f'{outdir}/04.Degenerate_Bases/primer2taxid.table', f'{outdir}/05.Primer_Reevaluate', info['nt_database'], info['database_path'])
            outshells.append(Primer_Reevaluate_shell)

        with open(f'{outdir}/00.shell/run{run_yaml}.sh', 'w') as SHELL1:
            run_shell = "\n".join(outshells)
            print(run_shell, file=SHELL1)

    if (first_step == '3'):
        outshells = []
        reference_dir = os.path.abspath(info['reference_dir'])
        chinese = info['chinese']
        conserved_bed = info['conserved_bed']

        Nested_PCR_shell = Nested_PCR(reference_dir, conserved_bed,
                                      f'{outdir}/03.Nested_PCR', info['name2taxid'], info['nt_database'], info['evalue'], chinese, info['base_num'])
        outshells.append(Nested_PCR_shell)

        if (int(last_step) >= 4):
            Degenerate_Bases_shell = Degenerate_Bases(
                f'{outdir}/03.Nested_PCR/Merge/{chinese}_merge_primer.xls', f'{reference_dir}/all.fna', f'{outdir}/04.Degenerate_Bases/')
            outshells.append(Degenerate_Bases_shell)

        if (int(last_step) >= 5):
            Primer_Reevaluate_shell = Primer_Reevaluate(
                f'{outdir}/04.Degenerate_Bases/primer2taxid.table', f'{outdir}/05.Primer_Reevaluate', info['nt_database'], info['database_path'])
            outshells.append(Primer_Reevaluate_shell)

        with open(f'{outdir}/00.shell/run{run_yaml}.sh', 'w') as SHELL1:
            run_shell = "\n".join(outshells)
            print(run_shell, file=SHELL1)

    if (first_step == '4'):
        outshells = []
        reference_dir = os.path.abspath(info['reference_dir'])
        chinese = info['chinese']
        merge_primer = info['merge_primer']

        Degenerate_Bases_shell = Degenerate_Bases(
            merge_primer, f'{reference_dir}/all.fna', f'{outdir}/04.Degenerate_Bases/')
        outshells.append(Degenerate_Bases_shell)

        if (int(last_step) >= 5):
            Primer_Reevaluate_shell = Primer_Reevaluate(
                f'{outdir}/04.Degenerate_Bases/primer2taxid.table', f'{outdir}/05.Primer_Reevaluate', info['nt_database'], info['database_path'])
            outshells.append(Primer_Reevaluate_shell)

        with open(f'{outdir}/00.shell/run{run_yaml}.sh', 'w') as SHELL1:
            run_shell = "\n".join(outshells)
            print(run_shell, file=SHELL1)

    if (first_step == '5'):
        outshells = []
        reference_dir = os.path.abspath(info['reference_dir'])
        primer2taxid_table = info["primer2taxid"]

        Primer_Reevaluate_shell = Primer_Reevaluate(
            primer2taxid_table, f'{outdir}/05.Primer_Reevaluate', info['nt_database'], info['database_path'])
        outshells.append(Primer_Reevaluate_shell)

        with open(f'{outdir}/00.shell/run{run_yaml}.sh', 'w') as SHELL1:
            run_shell = "\n".join(outshells)
            print(run_shell, file=SHELL1)

    if (first_step == '6'):
        pass


if __name__ == "__main__":
    Amplicon()
