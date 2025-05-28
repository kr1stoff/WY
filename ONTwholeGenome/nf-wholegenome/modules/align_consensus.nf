process BCFTOOLS_CONSENSUS {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.cons", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(vcf)
  tuple val(sample_id), path(ref_masked)

  output:
  tuple val(sample_id), path("consensus.fa")
  
  script:
  """
  bgzip -c $vcf > ${vcf}.gz
  bcftools index ${vcf}.gz
  cat $ref_masked | bcftools consensus ${vcf}.gz > consensus.fa
  """
}