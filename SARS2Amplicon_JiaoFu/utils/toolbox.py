#!/usr/bin/env python
import click
from pathlib import Path
import pandas as pd
import yaml
import re


@click.group(context_settings={'help_option_names': ['-h','--help']})
def cli():
    """[交付版本]二代测序新冠扩增子分析流程准备程序集."""

@cli.command()
@click.option('-x', '--excel', required=True, type=click.Path(exists=True), 
                help='信息输入表格. 模板见template/新冠样本信息表_生信.xlsx.')
@click.option('-d', '--outdir', required=True, help='输出结果目录.')
@click.option('-c', '--company', default='WY', show_default=True, type=click.Choice(['WY','KP']), help='哪个公司的报告模板.')
def get_input_yaml(excel, outdir, company):
    """获取input.yaml供主流程使用, 及sample_info.txt供WORD报告脚本使用."""
    Path(outdir).joinpath('materials').mkdir(exist_ok=True, parents=True)
    #input.yaml --> dict_yaml
    outdir_parent, outdir_name = Path(outdir).resolve().parent, Path(outdir).name
    dict_yaml = {
        'bed': '',
        'global_sars2': False,
        'library': outdir_name,
        'result_dir': str(outdir_parent),
        'company': company,
        'samples': {}}
    #sample_info.txt --> dict_sample_info
    df = pd.read_excel(excel)
    df.fillna('-', inplace=True)
    for row in df.iterrows():
        dict_sample_info = {
            'DanWei': row[1]['单位'],
            'KeShi': row[1]['科室'],
            'SampleNumber': row[1]['样本编号'],
            'ReportNumber': row[1]['报告编号'], 
            'SongJianDate': row[1]['送检日期'],
            'SampleType': row[1]['样本类型']
        }
        dict_yaml['samples'][row[1]['报告编号']] = [row[1]['FASTQ1']]
        if row[1]['FASTQ2'] != '-':
            dict_yaml['samples'][row[1]['报告编号']].append(row[1]['FASTQ2'])
        #output
        file_sample_info = f'{outdir}/materials/{row[1]["报告编号"]}.sample_info.txt'
        with open(file_sample_info, 'wt', encoding='utf-8', newline='') as g:
            for k, v in dict_sample_info.items():
                g.write(f'{k}\t{v}\n')
    #output
    with open(f'{outdir}/materials/input.yaml', 'wt', encoding='utf-8', newline='') as g:
        g.write(yaml.dump(dict_yaml, allow_unicode=True))

@cli.command()
@click.option('-d', '--outdir', required=True, help='新冠流程输出结果目录, 是输入也是输出目录.')
def get_quality_summary(outdir):
    """(SOLAR目录可用, 送检不行) 获取一个RUN的质控信息汇总表."""
    with open(f'{outdir}/cov_lineage/lineage_report_trans.xls') as f, \
        open(f'{outdir}/quality_summary.tsv', 'wt', encoding='utf-8', newline='') as g:
        g.write('\t'.join(['样本名', 'WHO命名', 'Pangolin谱系', 
        '比对率', '平均深度', '覆盖度', '深度≥10x', '深度≥30x', '深度≥100x', '均一性', 
        "过滤前总reads数", "过滤后总reads数", "过滤前总碱基数", "过滤后总碱基数", "过滤前序列平均长度", "过滤后序列平均长度", 
        "过滤前Q20", "过滤后Q20", "过滤前Q30", "过滤后Q30", "低质量reads数", "含N碱基过多的reads数", "低复杂度reads数", 
        "过滤前重复reads比例", "过滤前重复reads比例", "adapter"]) + '\n')
        next(f)
        for line in f:
            outlist = []
            llst = line.strip().split('\t')
            sample = llst[0]
            #pangolin 结果
            outlist.extend(llst[:3])
            #比对质量
            with open(f'{outdir}/{sample}/2.map/{sample}.bam_stats.txt') as d:
                next(d)
                outlist.extend(d.read().strip().split('\t'))
            #FQ质量
            with open(f'{outdir}/{sample}/1.qc/{sample}.basic.stat.txt') as d:
                next(d)
                for line in d:
                    if re.match('^[ATGCN]', line): continue #A碱基这种统计不要
                    values = line.strip().split('\t')[1:]
                    for val in values:
                        if val == '-': continue #过滤后没有"低质量reads数"之类的统计
                        outlist.append(val)
            g.write('\t'.join(outlist) + '\n')

if __name__ =='__main__':
    cli()
