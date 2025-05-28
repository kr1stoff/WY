import sys
from pathlib import Path
from subprocess import run
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
import pandas as pd
sys.path.append(Path(__file__).parent)
from .config import config


# 配置文件
params = config()


def faidx_get_fasta(rec, feat, fna, locus):
    """使用samtools_faidx提取locus区域序列, 写入SeqRecord."""
    faidx_location = f'{rec.id}:{feat.location.start+1}-{feat.location.end}' #221206 BCBio.GFF use python 0-base
    cml = f'{params.samtools} faidx {fna} {faidx_location}'
    res = run(cml, shell=True, capture_output=True, encoding='utf-8')
    contents = res.stdout.strip().split('\n')
    head = contents[0].replace('>','')
    seq = ''.join(contents[1:])
    orec = SeqRecord(Seq(seq), id=head, description=locus)
    return orec


def get_given_locus_tags(clustered_proteins, target_cluster):
    """获取蛋白簇的locus_tags"""
    with open(clustered_proteins) as f:
        for line in f:
            cluster, loca_strings = line.strip().split(':')
            loca = loca_strings.strip().split('\t')
            if cluster == target_cluster: break
    return loca


def write_clusters_fasta(loca, dir_fna, dir_gff, out_fa):
    """从各基因组中提取同源核心基因序列, 合并到一个FASTA文件."""
    records = []
    for locus in loca:
        prefix = '_'.join(locus.split('_')[:2])
        fna = list(Path(dir_fna).glob(f'{prefix}*.fna'))[0]
        gff = list(Path(dir_gff).glob(f'{prefix}*.gff'))[0]
        for rec in GFF.parse(open(gff)):
            if rec.features == []: continue
            for feat in rec.features:
                quali = feat.qualifiers
                if 'locus_tag' not in quali: continue
                locus_tag = quali['locus_tag'][0]
                if locus_tag == locus: break #找到目标
            if locus_tag == locus: break #找到目标
        orec = faidx_get_fasta(rec, feat, fna, locus)
        records.append(orec)
    SeqIO.write(records, out_fa, 'fasta')


def get_core_clusters(gene_presence_absence, thres_presence:float):
    """
    根据gene_presence_absence, 获取需要提取的核心蛋白簇列表
    gene_presence_absence   ::  roary结果文件gene_presence_absence.Rtab
    thres_presence          ::  出勤率阈值, 筛选出这个阈值以上的核心蛋白簇
    """
    df = pd.read_table(gene_presence_absence, sep='\t', encoding='utf-8', index_col=0)
    coregene_number_thres = int(df.shape[1]*thres_presence)
    rowsum = df.sum(1).sort_values(ascending=False)
    core_clusters = rowsum[rowsum>coregene_number_thres].index.values
    return core_clusters
