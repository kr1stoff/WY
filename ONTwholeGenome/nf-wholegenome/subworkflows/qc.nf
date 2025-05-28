include { NANOSTAT_PRE;NANOSTAT_POST;EXTRACT_NANOSTAT } from '../modules/nanostat.nf'
include { NANOFILT } from '../modules/nanofilt.nf'
include { NANOPLOT_PRE;NANOPLOT_POST } from '../modules/nanoplot.nf'


workflow QC {
  take:
  ch_fastq

  main:
  NANOSTAT_PRE(ch_fastq)
  NANOFILT(ch_fastq)
  NANOSTAT_POST(NANOFILT.out)
  NANOPLOT_PRE(ch_fastq)
  NANOPLOT_POST(NANOFILT.out)
  EXTRACT_NANOSTAT(NANOSTAT_PRE.out, NANOSTAT_POST.out)

  emit:
  fastq_cln = NANOFILT.out                // channel: [sample_id, clean_fastq]
}