process NANOFILT {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc", mode: "copy"
  conda params.conda_env.nano

  input:
  tuple val(sample_id), path(fastq)

  output: 
  tuple val(sample_id), path("${sample_id}.cln.fastq")

  script:
  """
  NanoFilt -q 10 --length 200 $fastq > ${sample_id}.cln.fastq
  """
}
