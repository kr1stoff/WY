process CANU {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(fastq)

  output: 
  tuple val(sample_id), path("canu/${sample_id}.contigs.fasta")

  script:
  """
  canu -p ${sample_id} -d canu genomeSize=20k maxThreads=${params.threads} -nanopore $fastq 
  """
}