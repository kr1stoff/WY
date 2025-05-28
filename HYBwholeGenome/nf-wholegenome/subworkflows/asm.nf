include { UNICYCLER } from '../modules/unicycler.nf'
include { QUAST } from '../modules/quast.nf'
include { CHECKM;CHECKM_PARSE } from '../modules/checkm.nf'
include { ASSEMBLY_DEPTH;GC_DEPTH_CONSISTENCY } from '../modules/asm_depth.nf'
include { ASSEMBLY_STAT } from '../modules/asm_stat.nf'
include { CHECKV } from '../modules/checkv.nf'


workflow ASSEMBLY {
  take:
  ch_ont_cln_fastq
  ch_ilmn_cln_fastq

  main:
  UNICYCLER(ch_ont_cln_fastq, ch_ilmn_cln_fastq)
  QUAST(UNICYCLER.out.fasta_asm)
  
  if ( params.biwf == "virassy") {
    CHECKV(UNICYCLER.out.fasta_asm)
  }
  else {
    CHECKM(UNICYCLER.out.dir_asm)
    CHECKM_PARSE(CHECKM.out)
  }
  
  ASSEMBLY_STAT(UNICYCLER.out.fasta_asm)
  ASSEMBLY_DEPTH(
    UNICYCLER.out.fasta_asm,
    ch_ont_cln_fastq,
    ch_ilmn_cln_fastq
  )
  GC_DEPTH_CONSISTENCY(ASSEMBLY_DEPTH.out.depth_asm, UNICYCLER.out.fasta_asm)

  emit:
  fasta_asm = UNICYCLER.out.fasta_asm            // channel: [sample_id, assembly_fasta]
  bam_asm = ASSEMBLY_DEPTH.out.bam_asm      // channel: [sample_id, assembly_bam]
}