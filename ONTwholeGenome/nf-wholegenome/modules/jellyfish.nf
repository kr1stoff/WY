process JELLYFISH_COUNT {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/jellyfish", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(fastq)

  output: 
  tuple val(sample_id), path("${sample_id}.jf")

  script:
  """
  jellyfish count -C -t ${params.threads} -m 21 -s 5G -o ${sample_id}.jf $fastq
  """
}


process JELLYFISH_HISTO {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/jellyfish", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(jf_count)

  output: 
  tuple val(sample_id), path("${sample_id}.histo")

  script:
  """
  jellyfish histo -t ${params.threads} -o ${sample_id}.histo $jf_count
  """
}
