process ISLAND {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred", mode: "copy"

  input:
  tuple val(sample_id), path(dir_predict)

  output:
  tuple val(sample_id), path("island", type:"dir")

  when:
  params.biwf == "batassy"
  
  script:
  """
  python ${projectDir}/bin/multi_Island.py -i ${dir_predict}/predict.gbk -o island
  """
}