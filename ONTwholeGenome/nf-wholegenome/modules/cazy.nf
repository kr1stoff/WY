process CAZY {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno", mode: "copy"
  conda params.conda_env.denovo

  input:
  tuple val(sample_id), path(dir_predict)
  
  output:
  tuple val(sample_id), path("CAZy", type:"dir")
  
  script:
  """
  run_dbcan ${dir_predict}/predict.faa protein --db_dir ${params.db.CAZY_dir}  --out_dir CAZy
  """
}


process CAZY_PARSE1 {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/CAZy", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(dir_cazy)

  output:
  tuple val(sample_id), path("CAZY.kind_stat.txt")
  path("CAZY.txt")
  
  script:
  """
  csvtk join -t -L -f 1 ${dir_cazy}/hmmer.out ${params.db.CAZY_dir}/CAZyDB.07302020.fam-activities-hmm.txt | \
    csvtk round -t -f 'Coverage' -n 2 | csvtk mutate -t -f 1 -p '^(\\D+)' -n 'kind' > CAZY.txt
  csvtk freq -t -f 'kind' CAZY.txt > CAZY.kind_stat.txt
  """
}


process CAZY_PARSE2 {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/CAZy", mode: "copy"

  input:
  tuple val(sample_id), path(stat_cazy)
  
  output:
  tuple val(sample_id), path("CAZy_anno_stats.png")
  
  script:
  """
  Rscript ${projectDir}/bin/draw_CAZy.R ${stat_cazy} CAZy_anno_stats.png
  """
}