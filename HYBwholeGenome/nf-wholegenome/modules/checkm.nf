process CHECKM {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.checkm

  input:
  tuple val(sample_id), path(dir_asm)
  
  output:
  tuple val(sample_id), path("checkm", type:"dir")
  
  script:
  """
  # FileNotFoundError: [Errno 2] No such file or directory: '~/.checkm/hmms/phylo.hmm'
  # https://github.com/Ecogenomics/CheckM/wiki/Installation#required-reference-data

  export CHECKM_DATA_PATH=${params.db.CheckM_dir}
  checkm lineage_wf -t ${params.threads} -q --tab_table -f checkm.txt -x fasta $dir_asm checkm
  """
}


process CHECKM_PARSE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/checkm", mode: "copy"

  input:
  tuple val(sample_id), path(dir_checkm)
  
  output:
  tuple val(sample_id), path("checkm_parse.txt")
  
  script:
  """
  python ${projectDir}/bin/parse_checkm_result.py -i ${dir_checkm}/storage/bin_stats_ext.tsv -o checkm_parse.txt
  """
}