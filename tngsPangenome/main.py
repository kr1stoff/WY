#!/usr/bin/env python
# @CreateTime       : 2022/12/05
# @Author           : mengxf
# @version          : v1.0.4
# @LastModified     : 2023/08/09

from pathlib import Path
import click
import logging
from modules.conserved_domain_from_roary import *
from modules.core_gene_fa import *
from modules.roary_stats import genome_coregene_count
from modules.batch_prokka import batch_prokka
from modules.accession_list import accession_list
from modules.predict_schedule import predict_schedule
from modules.degenerate import DegenerateMSA


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.group(context_settings={'help_option_names':['-h','--help']})
def cli():
    """WY扩增子平台引物设计泛基因组."""


@cli.command()
@click.option('-p', '--gene_presence_absence', required=True, type=click.Path(exists=True), 
              help='roary结果gene_presence_absence.Rtab文件.')
@click.option('-c', '--clustered_proteins', required=True, type=click.Path(exists=True), help='Roary结果clustered_proteins文件.')
@click.option('-g', '--ref_gff', required=True, type=click.Path(exists=True), help='参考/代表性基因组的.gff文件.')
@click.option('-o', '--out_table', default='core_cluster_location_presence.txt', show_default=True, 
              help='输出核心基因群坐标和出勤率表格.')
@click.option('-b', '--out_bed', default='core_cluster.bed', show_default=True, help='输出核心蛋白簇.bed文件.')
@click.option('-P', '--thres_presence', type=float, default=0.95, show_default=True, 
              help='出席率(presence)阈值, 在基因组集中的共有频率阈值.')
@click.option('-l', '--thres_length', type=int, default=100, show_default=True, help='locus_tag长度阈值, 小于该值的过滤掉.')
def get_conserved_from_roary(gene_presence_absence, clustered_proteins, ref_gff, out_table, out_bed, thres_presence, thres_length):
    """
    使用roary结果和gff文件坐标, 汇总表格.\n
    1. 获取核心蛋白群(core_cluster), 群名称及出勤率; 
    2. 在第1步核心蛋白群范围内, 获取cluster与locus对照表; 
    3. 在prokka注释结果中, 获取locus坐标信息; 
    4. 合并上面三部结果生成表格, cluster-locus-location-presence. 
    """
    #630 Yersinia enterocolitica 参考基因组assemble_accession: GCF_001160345.1
    #input & output
    ref_assembly = Path(ref_gff).stem
    #run
    dict_core_cluster_presence = get_core_cluster_presence_dict(gene_presence_absence, ref_assembly, thres_presence)
    dict_cluster_locus = get_cluster_locus_dict(clustered_proteins, dict_core_cluster_presence, ref_assembly)
    dict_locus_location = get_locus_locations_dict(ref_gff, thres_length)
    create_cluster_location_presence_table(
        out_table,
        out_bed, 
        dict_cluster_locus,
        dict_locus_location, 
        dict_core_cluster_presence
    )


@cli.command()
@click.option('-p', '--gene_presence_absence', required=True, type=click.Path(exists=True), 
              help='roary结果gene_presence_absence.Rtab文件.')
@click.option('-c', '--clustered_proteins', required=True, type=click.Path(exists=True), help='Roary结果clustered_proteins文件.')
@click.option('-f', '--dir_fna', required=True, type=click.Path(exists=True), help='物种下所有基因组FASTA文件目录.')
@click.option('-g', '--dir_gff', required=True, type=click.Path(exists=True), help='物种下所有基因组GFF文件目录.')
@click.option('-O', '--out_dir', default='multi_core_clusters', show_default=True, help='输出多基因组核心基因FASTA文件目录.')
@click.option('-P', '--thres_presence', type=float, default=0.95, show_default=True, 
              help='出席率(presence)阈值, 在基因组集中的共有频率阈值.')
def get_coregene_fa(gene_presence_absence, clustered_proteins, dir_fna, dir_gff, out_dir, thres_presence):
    """获取同物种下各基因组特定核心基因序列, 写入一个FASTA文件中. 后续多序列比对使用."""
    Path(out_dir).mkdir(exist_ok=True)
    for target_cluster in get_core_clusters(gene_presence_absence, thres_presence):
        loca = get_given_locus_tags(clustered_proteins, target_cluster)
        write_clusters_fasta(loca, dir_fna, dir_gff, Path(out_dir).joinpath(target_cluster+'.fna'))


@cli.command()
@click.option('-p', '--gene_presence_absence', required=True, type=click.Path(exists=True), 
              help='roary结果gene_presence_absence.Rtab文件.')
@click.option('-g', '--genome_coregene', default='genome_coregene_counts.txt', show_default=True, help='输出基因组核心基因数表.')
@click.option('-c', '--coregene_genome', default='coregene_genome_counts.txt', show_default=True, help='输出核心基因基因组数表.')
def get_extra_roary_stats(gene_presence_absence, genome_coregene, coregene_genome):
    """
    roary额外的统计\n
    1. gene_presence_absence.Rtab每个基因组包含核心基因簇的数量统计和每个核心基因簇的基因组数.
    """
    genome_coregene_count(gene_presence_absence, genome_coregene, coregene_genome)


@cli.command()
@click.option('-a', '--accessions', required=True, type=click.Path(exists=True), 
              help='accessions号文件, 子命令"get_accession_list"的输出结果.')
@click.option('-f', '--dir_fna', required=True, type=click.Path(exists=True), help='物种下所有基因组FASTA文件目录.')
@click.option('-g', '--dir_gff', required=True, type=click.Path(exists=True), help='物种下所有基因组GFF文件目录.')
def run_batch_prokka(dir_fna, dir_gff, accessions):
    """批量跑prokka, 输入为物种目录包含多个基因组FASTA, 输出为GFF物种目录."""
    batch_prokka(dir_fna, dir_gff, accessions)


@cli.command()
@click.option('-e', '--excel', required=True, type=click.Path(exists=True), help='钉钉项目进度表, 表头必有 分类/过滤完成/taxid.')
@click.option('-d', '--download_genome_dir', type=click.Path(exists=True), help='下载基因组的文件目录.')
@click.option('-w', '--work_dir', type=click.Path(exists=True), help='工作目录, prokka/roary/conserve都在这个目录做.')
@click.option('-k', '--skip_completed', is_flag=True, help='跳过已经生成accessions_file且不是空文件或已经完成预测的物种.')
def get_accession_list(excel, download_genome_dir, work_dir, skip_completed):
    """获取各物种accessions(已下载)列表,并软连接fna到'{taxid}/fnas'目录下."""
    accession_list(excel, download_genome_dir, work_dir, skip_completed)


@cli.command()
@click.option('-e', '--excel', required=True, type=click.Path(exists=True), help='钉钉项目进度表, 表头必有 分类/过滤完成/taxid.')
@click.option('-d', '--download_genome_dir', type=click.Path(exists=True), help='下载基因组的文件目录.')
@click.option('-w', '--work_dir', type=click.Path(exists=True), help='工作目录, prokka/roary/conserve都在这个目录做')
def get_predict_schedule(excel, download_genome_dir, work_dir):
    """获取预测进度表."""
    predict_schedule(excel, download_genome_dir, work_dir)


@cli.command()
@click.option('-p', '--primertable', required=True, type=click.Path(exists=True), help='引物表格文件,包含内左右引物/探针/外左右引物/参考登录号.')
@click.option('-a', '--allfna', required=True, type=click.Path(exists=True), help='单物种目录下的all.fna文件.')
@click.option('-d', '--outdir', default='degenerate', show_default=True, help='输出结果目录.')
def translate_degenerate_base(primertable, allfna, outdir):
    """
    !!!!!!!!! 已转移至 qpcrDesign 目录 !!!!!!!!! 
    
    简并碱基转换.
    """
    pipe = DegenerateMSA(primertable, allfna, outdir)
    pipe.degenerate()

if __name__ == '__main__':
    cli()
