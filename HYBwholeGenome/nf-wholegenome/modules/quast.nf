process QUAST {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(fasta_asm)
  
  output:
  tuple val(sample_id), path("quast", type:"dir")
  
  script:
  """
  quast --plots-format png -o quast -t ${params.threads} $fasta_asm
  """
}