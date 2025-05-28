#!/usr/bin/env bash

######## Options ########
printUsage() {
    echo -e 'PROGRAM <sample_name> <rawdata_path> <output_dir>'
    echo -e '\tsample_name    ::   Sample name. e.g. SRR12336750'
    echo -e '\trawdata_path   ::   Rawdata path, STRING separated by ",". e.g. PRJNA649101/SRR12336750_1.fastq.gz,PRJNA649101/SRR12336750_2.fastq.gz'
    echo -e '\toutput_dir     ::   Project output direcotory. e.g. "."\n'
    echo -e 'Example:  PROGRAM SRR12336750 PRJNA649101/SRR12336750_1.fastq.gz,PRJNA649101/SRR12336750_2.fastq.gz results/PRJNA649101'
    exit 1
}
#` check arguments number
if [ $# != 3 ];then
    echo -e 'Need 3 Arguments!\n'
    printUsage
fi

######## Path ########
#~ $1 sample_name, $2 rawdata_path, $3_output dir
export SAMPLE_NAME=$1
export RAWDATA=$2
export RESULT_DIR=$3
READ1=$(echo ${RAWDATA} | awk 'split($0, a, ","){print a[1]}')
READ2=$(echo ${RAWDATA} | awk 'split($0, a, ","){print a[2]}')
#~ check file exist
if test ! -f $READ1 -o ! -f $READ2;then
    echo -e 'No rawdata file!\n'
    printUsage
fi

echo "########## START ##########"
# VARIABLE
# RAWDATA, RESULT_DIR, SAMPLE_NAME
# THREADS, SCRIPTS_PATH

echo "########## INPUT ##########"
echo "SAMPLE_NAME: $SAMPLE_NAME"
echo "RAWDATA:     $RAWDATA"

######## Main ########
## file path tree
mkdir -p ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/fastqc
mkdir -p ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim
mkdir -p ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/map_qc
mkdir -p ${RESULT_DIR}/${SAMPLE_NAME}/3.variants
mkdir -p ${RESULT_DIR}/${SAMPLE_NAME}/4.consensus
mkdir -p ${RESULT_DIR}/${SAMPLE_NAME}/5.phylogenetic

# QC & Trim
echo "########## QC & Trim ##########"
#~ fastqc: 11.34s
fastqc -t ${THREADS} \
    --extract \
    ${READ1} ${READ2} \
    -o ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/fastqc
#~ fastp 
fastp -q 15 -u 40 -l 25 --thread $THREADS \
    --cut_right --cut_window_size 20 --cut_mean_quality 30 --correction \
    -i ${READ1} \
    -I ${READ2} \
    -o ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.clean.R1.fq \
    -O ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.clean.R2.fq \
    -j ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.json \
    -h ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.html
#` fastp json2table, output: fastq_stats.txt
python ${SCRIPTS_PATH}/utils/fastp_json2table.py ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.json

# Mapping
echo "########## Mapping ##########"
#~ bwa index: 4.5s
#bwa index ${REFERENCE}
#~ bwa mem: 58.1s
bwa mem -t $THREADS -M -R '@RG\tID:'$SAMPLE_NAME'\tSM:'$SAMPLE_NAME \
    ${REFERENCE} \
    ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.clean.R1.fq \
    ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.clean.R2.fq \
    | samtools view -@ 4 -hbS - \
    | samtools sort -@ 4 -o ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.sorted.bam - 
#~ samtools index: 1.9s 
samtools index ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.sorted.bam
#~ bam stats
#` samtools depth -a
samtools depth -a ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.sorted.bam \
    > ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/map_qc/${SAMPLE_NAME}.bam.depth
#` output dir same as input. 
# 1 bam_stats.txt, 2 genome_coverage_depth.png, 3 genome_coverage_depth_ylim1000.png
Rscript ${SCRIPTS_PATH}/utils/genome_coverage.R \
    ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/map_qc/${SAMPLE_NAME}.bam.depth
#~ [20220225 update] low coverage regions, region > 20bp
bedtools genomecov -bga \
    -ibam ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.sorted.bam \
    | awk -v cov="${MIN_COVERAGE}" '$4<cov' \
    | bedtools merge -i - \
    | awk '$3-$2>20' \
    1> ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.lowcovmask.bed

# Variants Calling
echo "########## Variants Calling ##########"
#~ maskfasta before consensus 
bedtools maskfasta \
    -fi ${REFERENCE} \
    -bed ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.lowcovmask.bed \
    -fo ${RESULT_DIR}/${SAMPLE_NAME}/4.consensus/${SAMPLE_NAME}.reference_masked.fa
#~ [20220215 update] freebayes min_cov > 10x
freebayes -F ${MIN_ALT_FRACTION} -C ${MIN_ALT_COUNT} \
    -f ${RESULT_DIR}/${SAMPLE_NAME}/4.consensus/${SAMPLE_NAME}.reference_masked.fa \
    ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.sorted.bam \
    > ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.vcf
#~ snpEff annotation, database NC_045512.2
snpEff NC_045512.2 \
    ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.vcf \
    -htmlStats ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/snpEff_summary.html \
    -csvStats ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/snpEff_summary.csv \
    > ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.snpeff.vcf
#` filter and translate. vcf2table output: trans.vcf
python ${SCRIPTS_PATH}/utils/vcf2table.py ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.snpeff.vcf
#~ snptree use snps.recode.vcf
#` vcftools: out SRR12336750.snps.recode.vcf
vcftools --vcf ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.vcf \
    --remove-indels --recode --recode-INFO-all \
    --out ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.snps
#` compress & bcftools index
bgzip -c ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.snps.recode.vcf \
    > ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.snps.recode.vcf.gz
bcftools index ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.snps.recode.vcf.gz

# Consensus Sequences Fasta
echo "########## Consensus ##########"
#~ filter mutations for nextclade (nextclade.vcf): 1remove indels
python ${SCRIPTS_PATH}/utils/vcf_filter_4nextclade.py ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.vcf
#~ index
bgzip -c ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.nextclade.vcf \
    > ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.nextclade.vcf.gz
bcftools index ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.nextclade.vcf.gz
# consensus
bcftools consensus -p ${SAMPLE_NAME}\  \
    --haplotype A \
    -f ${RESULT_DIR}/${SAMPLE_NAME}/4.consensus/${SAMPLE_NAME}.reference_masked.fa \
    -o ${RESULT_DIR}/${SAMPLE_NAME}/4.consensus/${SAMPLE_NAME}.consensus.fa \
    ${RESULT_DIR}/${SAMPLE_NAME}/3.variants/${SAMPLE_NAME}.nextclade.vcf.gz
sed -i 's/ .*//g' ${RESULT_DIR}/${SAMPLE_NAME}/4.consensus/${SAMPLE_NAME}.consensus.fa

# Single Sample Phylogenetic Tree
echo "########## Phylogenetic Tree ##########"
#~ merge consensus.fa
cat ${SC2ISOLATES} \
    ${RESULT_DIR}/${SAMPLE_NAME}/4.consensus/${SAMPLE_NAME}.consensus.fa \
    > ${RESULT_DIR}/${SAMPLE_NAME}/5.phylogenetic/${SAMPLE_NAME}.sc2isolates.fa
#~ multiple sequences alignment
mafft --thread 2 --auto --maxiterate 1000 \
    ${RESULT_DIR}/${SAMPLE_NAME}/5.phylogenetic/${SAMPLE_NAME}.sc2isolates.fa \
    > ${RESULT_DIR}/${SAMPLE_NAME}/5.phylogenetic/${SAMPLE_NAME}.aln.fa
#~ FastTree
fasttree -nt ${RESULT_DIR}/${SAMPLE_NAME}/5.phylogenetic/${SAMPLE_NAME}.aln.fa \
    > ${RESULT_DIR}/${SAMPLE_NAME}/5.phylogenetic/${SAMPLE_NAME}.tree
#~ multi-seq phylogenetic tree 
Rscript ${SCRIPTS_PATH}/utils/tree_for_SC2Pipe.R \
    ${RESULT_DIR}/${SAMPLE_NAME}/5.phylogenetic/${SAMPLE_NAME}.tree

# remove intermedia files
rm ${RESULT_DIR}/${SAMPLE_NAME}/1.qc/trim/${SAMPLE_NAME}.*.fq
# rm ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.sorted.bam*
# rm ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.trim_primer.bam
rm ${RESULT_DIR}/${SAMPLE_NAME}/2.mapping/${SAMPLE_NAME}.mpileup

echo "########## END ##########"