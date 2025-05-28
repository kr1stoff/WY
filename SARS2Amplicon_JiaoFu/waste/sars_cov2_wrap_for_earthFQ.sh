#!/usr/bin/env bash

##############################################################################
# Name:        SARS-Cov-2 Analysis Pipeline
# Version:     v1.0
# Author:      mengxf
#=============================================================================
# Example:
# bash scripts/sars_cov2_wrap.sh fq.input 2333
##############################################################################

printUsage() {
    echo -e 'PROGRAM <fq.input> <task id>'
    exit 1
}
if [ $# != 2 ];then
    echo -e 'Need 2 Argumengt!\n'
    printUsage
elif [ ! -f ${1} ];then
    echo -e 'No fq.input file!\n'
    printUsage
fi
# assign arguments
# FQINPUT=/sdbb/share/test_data/Sars_Cov2_Amplicon/fq.input
# TASKID=2333
FQINPUT=$1
TASKID=$2

######## Configure ########
# main scripts path (all relative path)
# export SCRIPTS_PATH="/sdbb/share/pipeline/Sars_Cov2_Amplicon"
export SCRIPTS_PATH=$(dirname $0)
#~ assignment arguments
source ${SCRIPTS_PATH}/libs/sars_cov2.env
#~ Activate Conda Environment
source ${CONDA_ACTIVATE} ${MAIN_CONDA_ENV}

## 1. intermedia files
#~ RESULT_DIR
export RESULT_DIR=${EARTH_ANALYSIS}/${TASKID}
python3 ${SCRIPTS_PATH}/utils/pre_process_for_earthFQ.py $FQINPUT $RESULT_DIR

## 2. Ctreate File Path
mkdir -p ${RESULT_DIR}/logs
mkdir -p ${RESULT_DIR}/intermedia
mkdir -p ${RESULT_DIR}/snp_phylogenetic_tree
mkdir -p ${RESULT_DIR}/multiseq_phylogenetic_tree
mkdir -p ${RESULT_DIR}/cov_lineage

## 3. parallel samples 
#parallel -N 2 -j $JOBS \
#    --xapply "bash ${SCRIPTS_PATH}/sars_cov2_pipe.sh {1} {2} ${RESULT_DIR} \
#                > ${RESULT_DIR}/logs/{1}_sc2pipe.o 2> ${RESULT_DIR}/logs/{1}_sc2pipe.e" \
#    ::: $(tr '\t' '\n' < ${RESULT_DIR}/intermedia/sample_path_table.txt)
python3 ${SCRIPTS_PATH}/utils/parallel.py \
    --sample_table ${RESULT_DIR}/intermedia/sample_path_table.txt \
    --result_dir ${RESULT_DIR} \
    > ${RESULT_DIR}/logs/parallel.o \
    2> ${RESULT_DIR}/logs/parallel.e

## 4. SNP phylogenetic tree
bash ${SCRIPTS_PATH}/utils/snp_tree.sh

## 5. Multi-sequences phylogenetic tree
bash ${SCRIPTS_PATH}/utils/multiseq_tree.sh

## 6. Cov Lineage
conda activate $PANGOLIN_ENV
$PANGOLIN ${RESULT_DIR}/multiseq_phylogenetic_tree/merged.consensus.fa \
    -o ${RESULT_DIR}/cov_lineage \
    > ${RESULT_DIR}/logs/pangolin.o \
    2> ${RESULT_DIR}/logs/pangolin.e
conda deactivate
python ${SCRIPTS_PATH}/utils/cov_lineage_list.py \
    ${RESULT_DIR}/cov_lineage/lineage_report.csv \
    > ${RESULT_DIR}/logs/cov_lineage.o \
    2> ${RESULT_DIR}/logs/cov_lineage.e
    
# 7. html report
python ${SCRIPTS_PATH}/utils/report.py $TASKID \
    > ${RESULT_DIR}/logs/report.o \
    2> ${RESULT_DIR}/logs/report.e

# Deactivate Conda Environment
conda deactivate
