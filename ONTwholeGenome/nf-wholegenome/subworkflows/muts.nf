include { LOFREQ;EXTRECT_LOFREQ } from '../modules/lofreq.nf'


workflow MUTATIONS {
  take:
  ch_bam
  ch_ref_masked

  main:
  LOFREQ(ch_bam, ch_ref_masked)
  EXTRECT_LOFREQ(LOFREQ.out)

  emit:
  vcf = LOFREQ.out                // channel: [sample_id, vcf]
}