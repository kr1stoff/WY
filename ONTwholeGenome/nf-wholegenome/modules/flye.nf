process FLYE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.asm

  input:
  val(genome_size)
  tuple val(sample_id), path(fastq)

  output: 
  tuple val(sample_id), path("flye", type:"dir"), emit: dir_asm
  tuple val(sample_id), path("flye/assembly.fasta"), emit: fasta_asm
  tuple val(sample_id), path("flye/assembly_info.txt"), emit: assembly_info

  script:
  """
  flye --nano-corr $fastq --out-dir flye --threads ${params.threads} --genome-size $genome_size --asm-coverage 200
  """
}