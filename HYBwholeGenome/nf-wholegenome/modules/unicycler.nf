process UNICYCLER {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(fastq)
  tuple val(sample_id), path(read1), path(read2)
  
  output:
  tuple val(sample_id), path("unicycler", type:"dir"), emit: dir_asm
  tuple val(sample_id), path("unicycler/assembly.fasta"), emit: fasta_asm

  script:
  """
  unicycler -t ${params.threads} -1 $read1 -2 $read2 -l $fastq -o unicycler
  """
}