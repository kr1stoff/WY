process MASKED_FASTA {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.align", mode: "copy"
  conda params.conda_env.basic
  
  input:
  tuple val(sample_id), path(bam)

  output:
  tuple val(sample_id), path("ref.masked.fa")
  
  script:
  """
  bedtools genomecov -bga -ibam $bam | awk '\$4<10' | bedtools merge -i - | awk '\$3-\$2>20' > low_depth_mask.bed
  bedtools maskfasta -fi $params.fileRefAbso -bed low_depth_mask.bed -fo ref.masked.fa
  """
}