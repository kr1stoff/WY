#!/usr/bin/env bash

##############################################################################
# Name:        SARS-Cov-2 Analysis Pipeline
# Version:     v1.0
# Author:      mengxf
#=============================================================================
# Example:
# bash scripts/sars_cov2_wrap.sh materials/in.json &> log/sars_cov2_wrap.log
##############################################################################

printUsage() {
    echo -e 'PROGRAM <IN.JSON>'
    exit 1
}
if [ $# != 1 ];then
    echo -e 'Just Need 1 Argumengt!\n'
    printUsage
elif [ ! -f ${1} ];then
    echo -e 'No IN.JSON file!\n'
    printUsage
fi
IN_JSON=$1
######## Configure ########
# main scripts path (all relative path)
#~scripts: include sars_cov2_pipe.sh, utils/vcf2phylip.py, utils/tree.R, sars_cov2.env ...
# export SCRIPTS_PATH="/sdbb/bioinfor/mengxf/TASKS/BeiJingJiKong_DongAo_20220128/scripts"
export SCRIPTS_PATH=$(dirname $0)
## 1. env setting & intermedia files
#~ python scripts: input.json, sample_path_table.txt
#~ [INPUT API] materials/in.json
WORKENV=$(python3 ${SCRIPTS_PATH}/utils/pre_process_for_FASTQ.py ${IN_JSON})
eval $WORKENV
#~ assignment arguments
source ${SCRIPTS_PATH}/libs/sars_cov2.env
## Activate Conda Environment
source ${CONDA_ACTIVATE} ${MAIN_CONDA_ENV}

## 2. Ctreate File Path
mkdir -p ${RESULT_DIR}/logs
mkdir -p ${RESULT_DIR}/snp_phylogenetic_tree
mkdir -p ${RESULT_DIR}/multiseq_phylogenetic_tree
mkdir -p ${RESULT_DIR}/cov_lineage

## 3. parallel samples 
parallel -N 2 -j $JOBS \
    --xapply "bash ${SCRIPTS_PATH}/sars_cov2_pipe.sh {1} {2} ${RESULT_DIR} \
                > ${RESULT_DIR}/logs/{1}_sc2pipe.o 2> ${RESULT_DIR}/logs/{1}_sc2pipe.e" \
    ::: $(tr '\t' '\n' < ${RESULT_DIR}/intermedia/sample_path_table.txt)

## 4. SNP phylogenetic tree
#~ merge VCFs
ls ${RESULT_DIR}/details/*/3.variants/*.snps.recode.vcf.gz \
    > ${RESULT_DIR}/intermedia/vcf_list.txt
bcftools merge -m snps -f PASS,. --force-samples --output-type v \
    --file-list ${RESULT_DIR}/intermedia/vcf_list.txt \
    -o ${RESULT_DIR}/snp_phylogenetic_tree/merged.vcf \
    > ${RESULT_DIR}/logs/bcltools_merge.o \
    2> ${RESULT_DIR}/logs/bcltools_merge.e
#~ vcf2phylip --> snp_phylogenetic_tree/merged.min4.phy
python ${SCRIPTS_PATH}/utils/vcf2phylip.py \
    -i ${RESULT_DIR}/snp_phylogenetic_tree/merged.vcf \
    --output-folder ${RESULT_DIR}/snp_phylogenetic_tree \
    > ${RESULT_DIR}/logs/vcf2phylip.o \
    2> ${RESULT_DIR}/logs/vcf2phylip.e
#~ FastTree
fasttree -nt ${RESULT_DIR}/snp_phylogenetic_tree/merged.min4.phy \
    > ${RESULT_DIR}/snp_phylogenetic_tree/merged.min4.tree \
    2> ${RESULT_DIR}/logs/snpfasttree.e
#~ snp phylogenetic tree    
Rscript ${SCRIPTS_PATH}/utils/tree.R \
    ${RESULT_DIR}/snp_phylogenetic_tree/merged.min4.tree \
    > ${RESULT_DIR}/logs/snptree.R.o \
    2> ${RESULT_DIR}/logs/snptree.R.e
    
## 5. Multi-sequences phylogenetic tree
#~ merge consensus .fa
cat $REFERENCE \
    ${RESULT_DIR}/intermedia/lookback.fa \
    ${RESULT_DIR}/details/*/4.consensus/*.consensus.fa \
    > ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.fa
#~ multiple sequences alignment
mafft --auto --maxiterate 1000 \
    ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.fa \
    > ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.aln.fa \
    2> ${RESULT_DIR}/logs/mafft.e
#~ FastTree
fasttree -nt ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.aln.fa \
    > ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.tree \
    2> ${RESULT_DIR}/logs/multiseqfasttree.e
#~ multi-seq phylogenetic tree 
Rscript ${SCRIPTS_PATH}/utils/tree.R \
    ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.tree \
    > ${RESULT_DIR}/logs/multiseqtree.R.o \
    2> ${RESULT_DIR}/logs/multiseqtree.R.e

## 6. Cov Lineage
$PANGOLIN ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.fa \
    -o ${RESULT_DIR}/cov_lineage \
    > ${RESULT_DIR}/logs/pangolin.o \
    2> ${RESULT_DIR}/logs/pangolin.e
python ${SCRIPTS_PATH}/utils/cov_lineage_list.py \
    ${RESULT_DIR}/cov_lineage/lineage_report.csv \
    > ${RESULT_DIR}/logs/cov_lineage.o \
    2> ${RESULT_DIR}/logs/cov_lineage.e
    
# 7. html report
python ${SCRIPTS_PATH}/utils/report.py ${IN_JSON} \
    > ${RESULT_DIR}/logs/report.o \
    2> ${RESULT_DIR}/logs/report.e
    
# Deactivate Conda Environment
conda deactivate

# remove waste
rm ${RESULT_DIR}/details/*/report.pl
