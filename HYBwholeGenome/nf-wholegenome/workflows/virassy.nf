include { QC } from '../subworkflows/qc.nf'
include { ASSEMBLY } from '../subworkflows/asm.nf'
include { PREDICTION } from '../subworkflows/pred.nf'
include { VIR_ANNOTATION } from '../subworkflows/anno.nf'


workflow VIRASSY {
  take:
  ch_fastq

  main:
  QC(ch_fastq)
  ASSEMBLY(QC.out.ont_fastq_cln, QC.out.ilmn_fastq_cln)
  PREDICTION(ASSEMBLY.out.fasta_asm, ASSEMBLY.out.bam_asm)
  VIR_ANNOTATION(PREDICTION.out.dir_pred, PREDICTION.out.depth_pred)
}