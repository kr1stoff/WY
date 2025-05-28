include { QC } from '../subworkflows/qc.nf'
include { ALIGNMENT } from '../subworkflows/align.nf'
include { MUTATIONS } from '../subworkflows/muts.nf'
include { CONSENSUS } from '../subworkflows/cons.nf'

workflow VIRMAP {
  take:
  ch_fastq

  main:
  QC(ch_fastq)
  ALIGNMENT(QC.out.fastq_cln)
  MUTATIONS(ALIGNMENT.out.bam, ALIGNMENT.out.fa)
  CONSENSUS(MUTATIONS.out.vcf, ALIGNMENT.out.fa)
}