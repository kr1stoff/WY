import pandas as pd


def genome_coregene_count(gene_presence_absence, genome_coregene, coregene_genome):
    """
    gene_presence_absence横向纵向和, 基因组包含核心基因数和核心基因包含基因组数.
    gene_presence_absence   ::  roary结果文件gene_presence_absence.Rtab
    genome_coregene         ::  输出基因组核心基因数表, df.sum(axis=0)
    coregene_genome         ::  输出核心基因基因组数表, df.sum(axis=0)
    """
    df = pd.read_table(gene_presence_absence, sep='\t', encoding='utf-8', index_col=0)
    colsum = df.sum(0).sort_values(ascending=False)
    rowsum = df.sum(1).sort_values(ascending=False)
    colsum.to_csv(genome_coregene, sep='\t', encoding='utf-8', header=False)
    rowsum.to_csv(coregene_genome, sep='\t', encoding='utf-8', header=False)
