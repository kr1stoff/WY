process CHECKV {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.checkm

  input:
  tuple val(sample_id), path(fasta_pol)

  output: 
  tuple val(sample_id), path("checkv/quality_summary.tsv")

  script:
  """
  checkv end_to_end -d ${params.db.CHECKV_db} -t ${params.threads} $fasta_pol checkv
  """
}