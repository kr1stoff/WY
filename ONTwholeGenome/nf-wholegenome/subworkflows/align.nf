include { MINIMAP2 } from '../modules/minimap2.nf'
include { MAP_DEPTH;COVERAGE } from '../modules/align_coverage.nf'
include { MASKED_FASTA } from '../modules/low_depth_mask.nf'


workflow ALIGNMENT {
  take:
  ch_fastq

  main:
  MINIMAP2(ch_fastq)
  MAP_DEPTH(MINIMAP2.out)
  COVERAGE(MAP_DEPTH.out)
  MASKED_FASTA(MINIMAP2.out)

  emit:
  bam = MINIMAP2.out              // channel: [sample_id, sorted_bam]
  fa = MASKED_FASTA.out           // channel: [sample_id, masked_fa]
}