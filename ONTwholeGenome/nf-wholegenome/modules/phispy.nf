process PHISPY {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred", mode: "copy"

  input:
  tuple val(sample_id), path(dir_predict)

  output:
  tuple val(sample_id), path("phispy", type:"dir")
  
  when:
  params.biwf == "batassy"

  script:
  """
  python ${projectDir}/bin/run_phispy.py -i ${dir_predict}/predict.gbk -o phispy
  """
}