process ASSEMBLY_STAT {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/stat", mode: "copy"

  input:
  tuple val(sample_id), path(fasta_asm)
  tuple val(sample_id), path(assembly_info)

  output:
  tuple val(sample_id), path("${sample_id}_fa.stat.txt"), path("${sample_id}.length.png")

  script:
  """
  python ${projectDir}/bin/ass_fa_stat.py -fa ${fasta_asm} -a ${assembly_info} -p ${sample_id}
  """
}