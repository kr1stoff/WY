process NANOPLOT_PRE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc/nanoplot", mode: "copy"

  input:
  tuple val(sample_id), path(read1), path(read2), path(fastq)

  output: 
  tuple val(sample_id), path("${sample_id}.raw.stat"), path("${sample_id}.raw.png")

  script:
  """
  python ${projectDir}/bin/nanopore_fastq_quality.py --input $fastq --stats_file ${sample_id}.raw.stat \
    --stats_chart ${sample_id}.raw.png
  """
}


process NANOPLOT_POST {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc/nanoplot", mode: "copy"

  input:
  tuple val(sample_id), path(fastq)

  output: 
  tuple val(sample_id), path("${sample_id}.cln.stat"), path("${sample_id}.cln.png")
  
  script:
  """
  python ${projectDir}/bin/nanopore_fastq_quality.py --input $fastq --stats_file ${sample_id}.cln.stat \
    --stats_chart ${sample_id}.cln.png
  """
}
