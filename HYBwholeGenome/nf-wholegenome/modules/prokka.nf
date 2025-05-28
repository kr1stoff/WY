process PROKKA {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(fasta_asm)
  
  output:
  tuple val(sample_id), path("predict", type:"dir")
  
  script:
  if (params.biwf == "virassy")
    kingdom = "Viruses"
  else
    kingdom = "Bacteria"
  """
  prokka $fasta_asm --prefix predict --cpus ${params.threads} --outdir predict --kingdom $kingdom --force --addgenes
  """
}


process PROKKA_FEATURE_CN {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred/predict", mode: "copy"

  input:
  tuple val(sample_id), path(dir_predict)

  output:
  tuple val(sample_id), path("predict_kind.txt")

  script:
  """
  python ${projectDir}/bin/chinese_prokka.py -i ${dir_predict}/predict.txt -o predict_kind.txt
  """
}