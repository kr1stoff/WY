process NANOFILT {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc/nanofilt", mode: "copy"
  conda params.conda_env.nano

  input:
  tuple val(sample_id), path(read1), path(read2), path(fastq)

  output: 
  tuple val(sample_id), path("${sample_id}.cln.long.fastq")

  script:
  """
  NanoFilt -q 10 --length 200 $fastq > ${sample_id}.cln.long.fastq
  """
}
