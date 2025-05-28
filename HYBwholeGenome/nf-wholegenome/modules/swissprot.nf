process SWISSPROT {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/swissprot", mode: "copy"

  input:
  tuple val(sample_id), path(dir_predict)
  
  output:
  tuple val(sample_id), path("swissprot_result.tsv")

  script:
  if (params.biwf == "virassy")
    db_select = "Viruses"
  else
    db_select = "Bacteria"
  """
  python ${projectDir}/bin/UniProt/swissprot_db.py analysing -t ${params.threads} -i ${dir_predict}/predict.faa \
    -d ${params.db.SWISSPROT_dir} --softw_align diamond --db_select ${db_select} -o swissprot_result.tsv
  """
}