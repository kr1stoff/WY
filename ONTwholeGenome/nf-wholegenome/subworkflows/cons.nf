include { BCFTOOLS_CONSENSUS } from '../modules/align_consensus.nf'


workflow CONSENSUS {
  take:
  ch_vcf
  ch_ref_masked

  main:
  BCFTOOLS_CONSENSUS(ch_vcf, ch_ref_masked)
}