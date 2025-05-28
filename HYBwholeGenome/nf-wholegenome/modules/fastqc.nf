process FASTQC_PRE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc/fastqc", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(read1), path(read2), path(fastq)

  output:
  tuple val(sample_id), path("before", type:"dir")

  script:
  """
  mkdir before
  fastqc -t 16 --extract -o before $read1 $read2
  """
}


process FASTQC_POST {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc/fastqc", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(read1), path(read2)

  output:
  tuple val(sample_id), path("after", type:"dir")
  
  script:
  """
  mkdir after
  fastqc -t 16 --extract -o after $read1 $read2 
  """
}