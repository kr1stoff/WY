process PRED_DEPTH {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred/predict", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(dir_predict)
  tuple val(sample_id), path(bam_asm)

  output:
  tuple val(sample_id), path("bed_bam.merge_depth_cov.txt")

  script:
  """
  grep -v '#' ${dir_predict}/predict.gff | grep 'CDS' | grep 'ID=' | cut -f 1,4,5,9 | awk -F';' '{print \$1}' | sed 's/ID=//g' > predict.bed
  bedtools coverage -mean -a predict.bed -b $bam_asm > bed_bam.merge.depth
  bedtools coverage -a predict.bed -b $bam_asm > bed_bam.merge.cov
  paste bed_bam.merge.depth bed_bam.merge.cov | cut -f 4,5,13 > bed_bam.merge_depth_cov.txt
  """
}


process PRED_LENGTH {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred/predict", mode: "copy"

  input:
  tuple val(sample_id), path(dir_predict)

  output:
  tuple val(sample_id), path("predict.length.png")

  script:
  """
  python ${projectDir}/bin/predict_fa_stat.py -fa ${dir_predict}/predict.ffn -p predict
  """
}