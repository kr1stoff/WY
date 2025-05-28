#!/usr/bin/env bash

#~ merge VCFs
ls ${RESULT_DIR}/*/3.variants/*.snps.recode.vcf.gz \
    > ${RESULT_DIR}/intermedia/vcf_list.txt
bcftools merge -m snps -f PASS,. --force-samples --output-type v \
    --file-list ${RESULT_DIR}/intermedia/vcf_list.txt \
    -o ${RESULT_DIR}/snp_phylogenetic_tree/merged.vcf \
    > ${RESULT_DIR}/logs/bcltools_merge.o \
    2> ${RESULT_DIR}/logs/bcltools_merge.e
#~ vcf2phylip --> snp_phylogenetic_tree/merged.min1.phy
python ${SCRIPTS_PATH}/utils/vcf2phylip.py -m 1 \
    -i ${RESULT_DIR}/snp_phylogenetic_tree/merged.vcf \
    --output-folder ${RESULT_DIR}/snp_phylogenetic_tree \
    > ${RESULT_DIR}/logs/vcf2phylip.o \
    2> ${RESULT_DIR}/logs/vcf2phylip.e
#~ FastTree
fasttree -nt ${RESULT_DIR}/snp_phylogenetic_tree/merged.min1.phy \
    > ${RESULT_DIR}/snp_phylogenetic_tree/merged.min1.tree \
    2> ${RESULT_DIR}/logs/snpfasttree.e
#~ snp phylogenetic tree    
Rscript ${SCRIPTS_PATH}/utils/tree.R \
    ${RESULT_DIR}/snp_phylogenetic_tree/merged.min1.tree \
    > ${RESULT_DIR}/logs/snptree.R.o \
    2> ${RESULT_DIR}/logs/snptree.R.e
