process RGI {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/card", mode: "copy"
  conda params.conda_env.rgi

  input:
  tuple val(sample_id), path(dir_predict)
  
  output:
  tuple val(sample_id), path("card.txt")

  script:
  """
  rgi main -n ${params.threads} --input_type contig --clean -i ${dir_predict}/predict.ffn --output_file card
  """
}


process RGI_PARSE1 {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/card", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(dir_predict)
  tuple val(sample_id), path(card_txt)

  output:
  tuple val(sample_id), path("cut_id_card.txt")
  path("card.fa")
  
  script:
  """
  cut -f 1 card.txt | cut -f1 -d ' '| sed 's/ORF_ID/Contig/g' | sed 's/_1//g' > contig.id
  cut -f 8,10,11,15,16,21 card.txt > cut_card.txt
  paste contig.id cut_card.txt | csvtk join -tf 'ARO' - ${params.db.CARD_tab} > cut_id_card.txt
  seqkit grep -nrf contig.id ${dir_predict}/predict.ffn > card.fa
  """
}


process RGI_PARSE2 {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/card", mode: "copy"

  input:
  tuple val(sample_id), path(depth_predict)
  tuple val(sample_id), path(cut_id_card)

  output:
  tuple val(sample_id), path("detail_card_resistance.txt"), path("detail_card_resistance.xlsx")

  script:
  """
  python ${projectDir}/bin/parse.rgi.py $cut_id_card $depth_predict detail_card_resistance
  """
}