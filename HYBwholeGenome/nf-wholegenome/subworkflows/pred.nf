include { PROKKA;PROKKA_FEATURE_CN } from '../modules/prokka.nf'
include { PRED_DEPTH;PRED_LENGTH } from '../modules/pred_depth.nf'
include { PHISPY } from '../modules/phispy.nf'
include { ISLAND } from '../modules/island.nf'
include { REPEATMASKER;REPEATMASKER_PARSE } from '../modules/repeatmasker.nf'
include { TRF;TRF_PARSE } from '../modules/trf.nf'


workflow PREDICTION {
  take:
  ch_asm_fasta
  ch_asm_bam

  main:
  PROKKA(ch_asm_fasta)
  PROKKA_FEATURE_CN(PROKKA.out)
  PRED_DEPTH(PROKKA.out, ch_asm_bam)
  PRED_LENGTH(PROKKA.out)
  PHISPY(PROKKA.out)
  ISLAND(PROKKA.out)
  REPEATMASKER(ch_asm_fasta) | REPEATMASKER_PARSE
  TRF(ch_asm_fasta) | TRF_PARSE

  emit:
  dir_pred = PROKKA.out               // Return Channel: sample_id, dir_predict
  depth_pred = PRED_DEPTH.out
}