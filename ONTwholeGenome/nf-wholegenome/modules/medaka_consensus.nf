// Nanopore assembly genome polish
process MEDAKA_CONSENSUS {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.medaka

  input:
  tuple val(sample_id), path(fastq)
  tuple val(sample_id), path(fasta_asm)

  output: 
  tuple val(sample_id), path("medaka_consensus/consensus.fasta")

  script:
  """
  medaka_consensus -i $fastq -d $fasta_asm -o medaka_consensus -t ${params.threads}
  """
}