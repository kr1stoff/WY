#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/16 16:47
# @Last Modified by:   Ming
# @Last Modified time: 2022/6/16 16:47
import argparse
import datetime as datetime
import os
import subprocess
import sys

import dendropy
import pandas as pd
from Bio import SeqIO

args = argparse.ArgumentParser(description='Noro Virus Subtyping')
args.add_argument('--process', '-p', required=True, type=str, help='Options: makedb/subtype')
args.add_argument('--input', '-i', required=False, type=str, help='Input your sequence in FASTA format')
args.add_argument('--database', '-db', required=True, type=str, help='Genbank file of your database')
args.add_argument('--mapping', '-m', required=False, type=str, help='mapping genotype/serotype to accession id')
args.add_argument('--threshold', '-t', required=False, type=float, default=0.7,
                  help='Threshold of bootstrap value for assignment')
args.add_argument('--bootstrap', '-b', required=False, type=int, default=1000, help='Bootstrap times')
args.add_argument('--minimal', '-min', required=False, type=int, default=200, help='Minimal query length')
args.add_argument('--genes', '-g', required=False, type=str, nargs='*', help='Gene/CDS region you want to blast')
args.add_argument('--category', '-c', required=False, type=str, default=None, help='category to divide fasta/database')
args.add_argument('--source', '-s', required=False, type=str, default='all',
                  help='Assgin a smaller DB to blast taxonomy')
args.add_argument('--ambigouous', '-a', required=False, action='store_true', help='Ambigouously matching gene names')
args.add_argument('--separate', '-sp', required=False, action='store_true', help='Building seperate DB by category')
# fasta input??

args = args.parse_args()
print(
    f'Choosing process {args.process}, using database {args.database}\n-----------------------------------------------------------')
print(f'Using genes {args.genes}')


# sys.exit(1)
def makedb(gbk=args.database, mapping=args.mapping, categ_col=args.category, genes=args.genes, source=args.source):
    if mapping == None: print('Mapping file invalid'); sys.exit(1)
    mapping_geno = pd.read_csv(mapping, sep='\t')
    grep_cmd = fr"""less {gbk} |grep -E 'product="|gene="'|cut -d '"' -f2 | sed 's/^/"/'| sed 's/$/"/'| sort |uniq -c|sed 's/ "/@"/'| tr -d ' '|sed -r 's/(.*)@(.*)/\2\t\1/' > {gbk}-gene.list"""
    subprocess.run(grep_cmd, shell=True)
    print(f'Available gene regions for SUBTYPE are in file "{gbk}-gene.list"')
    gene_match_list = {}
    if genes is not None:
        if args.ambigouous == True:
            avail_gene_list = list(pd.read_csv(f'{gbk}-gene.list', header=None, sep='\t').loc[:, 0])
            for gene_chr in genes:
                gene_chr_list = []
                for gene in gene_chr.split(
                        '#'):  ###! support input form like :  "RdRp#Polymerase#dependent" "VP1#VP2#capsid"
                    gene_chr_list = gene_chr_list + list(filter(lambda v: gene.lower() in v.lower(), avail_gene_list))
                    gene_chr_list = list(set(gene_chr_list))
                gene_match_list[gene_chr] = gene_chr_list
        else:
            for gene_chr in genes:
                gene_match_list[gene_chr] = gene_chr.split('#')

    def gbk2fasta(gene=None):
        recs = SeqIO.parse(gbk, format='gb')
        dash_gene = f'-{gene}'.replace(' ', '_') if gene is not None else ''
        output_fasta = f'{gbk}-{source}{dash_gene}.fasta'
        if os.path.isfile(output_fasta): os.remove(output_fasta)
        acession_list = list(map(lambda x: x.strip(), list(mapping_geno.iloc[:, 0])))
        acession_list = [str(x) for x in acession_list]
        with open(output_fasta, 'a') as f:
            for rec in recs:
                mapping_pos = acession_list.index(rec.id.split('.')[0])
                genotype = list(mapping_geno.iloc[:, 1])[mapping_pos]
                categ = list(mapping_geno.loc[:, categ_col])[mapping_pos] if categ_col in list(
                    mapping_geno.columns) else 'Not_available'
                if str(categ) == 'nan': categ = 'Not_available'
                if str(genotype) == 'nan': genotype = 'Genotype_Null'
                anno_date = str(datetime.datetime.strptime(rec.annotations['date'], '%d-%b-%Y').date())
                fasta_id = rec.id.split('.')[
                               0] + '|' + anno_date + '|' + genotype + '|' + categ  # add geno sero
                if gene is not None:
                    for feat in rec.features:
                        if (feat.type == 'mat_peptide' or feat.type == 'CDS') and \
                                ((('product' in [key for key in feat.qualifiers]) and feat.qualifiers['product'][0] in
                                  gene_match_list[gene]) or \
                                 (('gene' in [key for key in feat.qualifiers]) and feat.qualifiers['gene'][0] in
                                  gene_match_list[gene])):
                            if feat.strand == 1:
                                fasta_seq = rec.seq[feat.location.start:feat.location.end]
                            if feat.strand != 1:
                                fasta_seq = rec.seq.reverse_complement()[feat.location.start:feat.location.end]
                            f.writelines(
                                '>' + feat.qualifiers['product'][0].replace(' ', '-') + '@' + fasta_id + '\n' + str(
                                    fasta_seq) + '\n')
                else:
                    f.writelines('>' + fasta_id + '\n' + str(rec.seq) + '\n')
            f.close()
        ref_align_cmd = f'mafft --auto --thread -1 {output_fasta} > {output_fasta}.fa'
        subprocess.run(ref_align_cmd, shell=True)

    def buildDB_by_categ(gbk, categ, gene=None):
        dash_gene = f'-{gene}'.replace(' ', '_') if gene is not None else ''
        dash_categ = f'-{categ}'.replace(' ', '_') if categ is not None else ''

        if categ is not None:
            if args.separate == True:
                with open(f'{gbk}{dash_categ}{dash_gene}.fasta', 'w') as f:
                    recs = SeqIO.parse(f'{gbk}-{source}{dash_gene}.fasta', format='fasta')
                    for rec in recs:
                        if categ in rec.id:
                            f.writelines('>' + rec.id + '\n' + str(rec.seq) + '\n')
                f.close()
                if os.stat(f'{gbk}{dash_categ}{dash_gene}.fasta').st_size == 0: return (None)
                cmd_BuildDB_categ = f'makeblastdb -in {gbk}{dash_categ}{dash_gene}.fasta -dbtype nucl -title {gbk}{dash_categ}{dash_gene} -out {gbk}{dash_categ}{dash_gene}.db'
                subprocess.run(cmd_BuildDB_categ, shell=True)
        else:
            cmd_BuildDB_categ = f'makeblastdb -in {gbk}-{source}{dash_gene}.fasta -dbtype nucl -title {gbk}-{source}{dash_gene} -out {gbk}-{source}{dash_gene}.db'
            if os.stat(f'{gbk}-{source}{dash_gene}.fasta').st_size == 0: return (None)
            print(f'---\nbuiding database for {gbk}-{source}{dash_gene}')
            subprocess.run(cmd_BuildDB_categ, shell=True)

    gbk2fasta(gene=None)
    buildDB_by_categ(gbk, None, None)

    if genes:
        for gene in genes:
            gbk2fasta(gene=gene)
            buildDB_by_categ(gbk, None, gene)

    if categ_col in list(mapping_geno.columns):
        categ_set = set(mapping_geno.loc[:, categ_col])
        for categ in categ_set:
            if str(categ) != 'nan':
                buildDB_by_categ(gbk, categ)
                if genes:
                    for gene in genes:
                        buildDB_by_categ(gbk, categ, gene)


def subtype(gbk=args.database, query=args.input, gene=None, categ_col=args.category, min_len=args.minimal,
            threshold=args.threshold, source='all'):
    def concatenate(query=args.input):
        query_list = list(SeqIO.parse(query, "fasta"))
        query_id = query_list[0].name
        if len(query_list) > 1:
            query_seq = ''.join([str(rec.seq) for rec in query_list])
            with open('concat.fasta', 'w') as f:
                f.writelines('>' + query_id + '\n' + query_seq)
                query = 'concat.fasta'
        else:
            query_seq = query_list[0].seq
        return query_id, query_seq, query

    query_id, query_seq, query = concatenate()
    dash_gene = f'-{gene}'.replace(' ', '_') if gene is not None else ''
    database = f'{gbk}-{source}{dash_gene}.db'

    # if query too short
    if len(query_seq) < min_len:
        print("Input sequence too short, results may not be reliable!")
        return (None)

    def taxonomy(gbk=gbk, query=query, database=database):
        cmd_blastn = f"blastn -query {query} -db {database} -out query.blast.res -outfmt 6"
        if os.path.exists(f'{database}.nhr') == False: print(f'{database} not found!!!');makedb(mapping=args.mapping)
        subprocess.run(cmd_blastn, shell=True)

        # if blast result empty
        if os.stat("query.blast.res").st_size == 0:
            print(f'Not a "{gbk}" "{gene}" sequence, genotype/serotype not assigned!')
            return (None)
        blast_res = pd.read_csv('query.blast.res', sep="\t", header=None)
        # if e-value too big
        min_e_value = list(blast_res[10]).index(min(blast_res[10]))
        if min_e_value > 1E-5: return (None)
        best_accession = \
            blast_res.loc[[list(blast_res[10]).index(min(blast_res[10]))], :][1][0].split('|')[0].split('@')[-1]
        best_taxon = ['/'.join(rec.annotations['taxonomy']) for rec in SeqIO.parse(gbk, "gb") if
                      rec.id.split('.')[0] == best_accession][0]
        print(f'---\nBest matching reference is {best_accession} \n')
        print('With taxonomy: ' + best_taxon)
        result_list.append(best_accession)
        result_list.append(best_taxon)
        return best_accession

    best_accession = taxonomy()
    if best_accession is None: return (None)

    def find_neighbour(ref_aligned='merge.fasta', query=query, query_id=query_id, threshold=threshold,
                       times=args.bootstrap):
        print(query, ref_aligned)
        cmd_align = f"mafft --auto --thread -1 --add {query} {ref_aligned} > merge.fa"
        subprocess.run(cmd_align, shell=True, stderr=subprocess.DEVNULL).stdout
        rec_list = list(SeqIO.parse('merge.fa', 'fasta'))
        if len(rec_list) <= 2: print(
            '\nToo less data in fasta file! Retry after change another "gene" or "database"!\n'); return (None)
        query = rec_list[0]
        query_start = len(query.seq) - len(query.seq.lstrip('-'))
        query_end = len(query.seq.rstrip('-'))
        with open('merge-truc.fa', 'w') as f:  # truncate
            for rec in rec_list:
                fasta_id = rec.id
                fasta_seq = rec.seq[query_start:query_end]
                f.writelines('>' + fasta_id + '\n' + fasta_seq + '\n')
            f.close()
        cmd_tree = f"rapidnj -n -t d -i fa merge-truc.fa -b {times} -x {ref_aligned}.nwk"
        print('\nBuilding tree...')
        subprocess.run(cmd_tree, shell=True).stdout

        print(f'{ref_aligned}.nwk')
        nwk = dendropy.Tree.get(path=f'{ref_aligned}.nwk', schema='newick')

        query_node = nwk.find_node_with_taxon_label(query_id)
        query_node_parent = query_node.parent_node
        support = float(query_node_parent.label) / 100 if query_node_parent.label is not None else 0
        result_list.append(support)
        sub_nwk = dendropy.Tree(seed_node=query_node_parent)

        if support >= threshold:

            print('---\nBootstrap support for the query node is ' + str(support) + '.')
            query_pos = sub_nwk.nodes().index(query_node)
            mtrx = sub_nwk.node_distance_matrix()
            dist_list = list(map(lambda v: mtrx.distance(sub_nwk.nodes()[query_pos], v), sub_nwk.nodes()))
            dist_list[query_pos] = max(dist_list) + 1
            dist_list = [(max(dist_list) + 1) if sub_nwk.nodes()[x].taxon == None else dist_list[x] for x in
                         range(0, len(dist_list))]  # convert non-tip length to max + 1

            import numpy as np
            nb_list = list(list(np.where(np.array(dist_list) == min(dist_list)))[0])
            nb_names = [str(sub_nwk.nodes()[x].taxon).replace("'", "") for x in nb_list if
                        sub_nwk.nodes()[x].taxon != None]
            print('---\nBest matching reference/genotype/serotype on the tree: \n')
            [print(x) for x in nb_names if not None] if len(nb_names) > 1 else print(nb_names[0])
            genotype = [i.split('|')[2] for i in nb_names]
            categ = [i.split('|')[-1] for i in nb_names]
            if len(nb_names) == 1: nb_names = nb_names[0]
            if len(set(genotype)) == 1: genotype = genotype[0]
            if len(set(categ)) == 1: categ = categ[0]
            result_list.append(nb_names)
            result_list.append(genotype)
            if categ_col != None:
                result_list.append(categ)
            else:
                result_list.append('-')

        else:
            print('---\nBootstrap support for the query node is ' + str(
                support) + f', below bootstrap threshold {threshold}.')

    def genotyping(gbk=gbk, best_accession=best_accession):
        categ = '-'
        if gene is not None: print(f'\nUsing gene "{gene}" ({gbk}-{source}.fasta.fa) to build tree...\n')
        if categ_col is not None:
            recs = SeqIO.parse(f'{gbk}-{source}.fasta', format='fasta')
            best_ref = [rec.id for rec in recs if best_accession in rec.id][0]
            categ = best_ref.split('|')[-1]  # categ might be unavailable
            print(f'---\n{categ_col} is {categ}, according to reference: {best_ref}\n')

        if args.separate == True and categ_col is not None:
            if os.path.exists(f'{gbk}-{categ}{dash_gene}.fasta') == False or os.stat(
                    f'{gbk}-{categ}{dash_gene}.fasta').st_size == 0:
                print('\nBuild blast DB for each "category" first using "makedb" !!!\n')
                sys.exit(1)
            find_neighbour(ref_aligned=f'{gbk}-{categ}{dash_gene}.fasta.fa')
        else:
            find_neighbour(ref_aligned=f'{gbk}-{source}{dash_gene}.fasta.fa')

    genotyping()


if args.process == 'makedb': makedb()

if args.process == 'subtype':
    genes = args.genes
    if not genes: genes = '-'
    with open('result.tsv', 'w') as OUT:
        header = ["Gene_name", "Genotype", args.category]
        print(*header, sep="\t", file=OUT)
        flag = 0
        for i in range(len(genes)):
            print(f'--------------------------------------------\n---\n---\nStarting gene "{genes[i]}" \n---\n')
            result_list = []
            if genes[i] != '-':
                subtype(gene=genes[i])
            else:
                subtype()
            while len(result_list) < 6: result_list.append('-')
            result_list.append(args.bootstrap)
            gene_name = genes[i]
            if flag == 0:
                result_list[5] = "NA"
            if flag == 1:
                result_list[4] = "NA"
            content = [gene_name, result_list[4], result_list[5]]
            print(*content, sep="\t", file=OUT)
            flag += 1
