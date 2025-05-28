include { QC } from '../subworkflows/qc.nf'
include { VIRUS_ASSEMBLY } from '../subworkflows/asm.nf'
include { PREDICTION } from '../subworkflows/pred.nf'
include { VIR_ANNOTATION } from '../subworkflows/anno.nf'


workflow VIRASSY {
  take:
  ch_fastq

  main:
  QC(ch_fastq)
  VIRUS_ASSEMBLY(QC.out.fastq_cln)
  PREDICTION(VIRUS_ASSEMBLY.out.fasta_asm, VIRUS_ASSEMBLY.out.bam_asm)
  VIR_ANNOTATION(PREDICTION.out.dir_pred, PREDICTION.out.depth_pred)
}