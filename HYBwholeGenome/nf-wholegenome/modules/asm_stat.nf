process ASSEMBLY_STAT {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/stat", mode: "copy"

  input:
  tuple val(sample_id), path(fasta_asm)

  output:
  tuple val(sample_id), path("${sample_id}_fa.stat.txt"), path("${sample_id}.length.png")
  
  script:
  """
  python ${projectDir}/bin/ass_fa_stat.py -fa ${fasta_asm} -p ${sample_id}
  """
}