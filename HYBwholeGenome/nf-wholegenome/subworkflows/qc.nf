include { NANOSTAT_PRE;NANOSTAT_POST;EXTRACT_NANOSTAT } from '../modules/nanostat.nf'
include { NANOFILT } from '../modules/nanofilt.nf'
include { NANOPLOT_PRE;NANOPLOT_POST } from '../modules/nanoplot.nf'
include { FASTQC_PRE;FASTQC_POST } from '../modules/fastqc.nf'
include { FASTP;EXTRACT_FASTP } from '../modules/fastp.nf'


workflow QC {
  take:
  ch_fastq

  main:
  // Nanopore
  NANOSTAT_PRE(ch_fastq)
  NANOPLOT_PRE(ch_fastq)
  NANOFILT(ch_fastq)
  NANOSTAT_POST(NANOFILT.out)
  NANOPLOT_POST(NANOFILT.out)
  EXTRACT_NANOSTAT(NANOSTAT_PRE.out, NANOSTAT_POST.out)
  // Illumina
  FASTQC_PRE(ch_fastq)
  FASTP(ch_fastq)
  EXTRACT_FASTP(FASTP.out.fastp_stat)
  FASTQC_POST(FASTP.out.clean_fastq)

  emit:
  ont_fastq_cln = NANOFILT.out                // channel: [sample_id, clean_fastq]
  ilmn_fastq_cln = FASTP.out.clean_fastq      // channel: [sample_id, clean_read1, clean_read2]
}