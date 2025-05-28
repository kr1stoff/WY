process MINIMAP2 {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.align", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(fastq)
  
  output:
  tuple val(sample_id), path("sort.bam")
  
  script:
  """
  minimap2 -d ${params.fileRefAbso}.mmi $params.fileRefAbso
  samtools faidx $params.fileRefAbso
  minimap2 -a -x map-ont -t ${params.threads} ${params.fileRefAbso}.mmi $fastq \
    | samtools view -@ ${params.threads} -bS -F 4 - | samtools sort -@ ${params.threads} -o sort.bam - 
  """
}
