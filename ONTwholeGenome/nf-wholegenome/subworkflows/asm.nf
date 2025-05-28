include { JELLYFISH_COUNT;JELLYFISH_HISTO } from '../modules/jellyfish.nf'
include { GENOMESCOPE2;GENOMESIZE1;GENOMESIZE2 } from '../modules/genomescope2.nf'
include { FLYE } from '../modules/flye.nf'
include { QUAST } from '../modules/quast.nf'
include { CHECKM;CHECKM_PARSE } from '../modules/checkm.nf'
include { ASSEMBLY_DEPTH;GC_DEPTH_CONSISTENCY } from '../modules/asm_depth.nf'
include { ASSEMBLY_STAT } from '../modules/asm_stat.nf'


workflow ASSEMBLY {
  take:
  ch_fastq_cln

  main:
  // 基因组大小预测
  JELLYFISH_COUNT(ch_fastq_cln)
  JELLYFISH_HISTO(JELLYFISH_COUNT.out)
  GENOMESCOPE2(JELLYFISH_HISTO.out)
  GENOMESIZE1(GENOMESCOPE2.out)
  GENOMESIZE2(GENOMESIZE1.out)
  // 组装及评估
  FLYE(GENOMESIZE2.out, ch_fastq_cln)
  QUAST(FLYE.out.fasta_asm)
  CHECKM(FLYE.out.dir_asm)
  CHECKM_PARSE(CHECKM.out)
  ASSEMBLY_DEPTH(FLYE.out.fasta_asm, ch_fastq_cln)
  GC_DEPTH_CONSISTENCY(ASSEMBLY_DEPTH.out.depth_asm, FLYE.out.fasta_asm)
  ASSEMBLY_STAT(FLYE.out.fasta_asm, FLYE.out.assembly_info)

  emit:
  fasta_asm = FLYE.out.fasta_asm            // channel: [sample_id, assembly_fasta]
  bam_asm = ASSEMBLY_DEPTH.out.bam_asm      // channel: [sample_id, assembly_bam]
}


include { CANU } from '../modules/canu.nf'
include { MEDAKA_CONSENSUS } from '../modules/medaka_consensus.nf'
include { CHECKV } from '../modules/checkv.nf'


workflow VIRUS_ASSEMBLY {
  take:
  ch_fastq_cln

  main:
  CANU(ch_fastq_cln)
  MEDAKA_CONSENSUS(ch_fastq_cln, CANU.out)
  CHECKV(MEDAKA_CONSENSUS.out)
  QUAST(CANU.out)
  ASSEMBLY_DEPTH(CANU.out, ch_fastq_cln)
  GC_DEPTH_CONSISTENCY(ASSEMBLY_DEPTH.out.depth_asm, CANU.out)
  ASSEMBLY_STAT(CANU.out)

  emit:
  fasta_asm = CANU.out                      // channel: [sample_id, assembly_fasta]
  bam_asm = ASSEMBLY_DEPTH.out.bam_asm      // channel: [sample_id, assembly_bam]
}