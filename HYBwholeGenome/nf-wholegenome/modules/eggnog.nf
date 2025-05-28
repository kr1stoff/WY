process EMAPPER {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno", mode: "copy"
  conda params.conda_env.denovo

  input:
  tuple val(sample_id), path(dir_predict)
  
  output:
  tuple val(sample_id), path("eggmapper", type:"dir")

  // emapper.py 需要提前创建输出目录
  script:
  """
  mkdir eggmapper
  emapper.py -m diamond --cpu ${params.threads} -d bact --override -i ${dir_predict}/predict.faa \
    --output $sample_id --output_dir eggmapper --data_dir ${params.db.EGGNOG_dir}
  cp -f ${dir_predict}/predict.faa eggmapper
  """
}


process EMAPPER_PARSE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno", mode: "copy"

  input:
  tuple val(sample_id), path(dir_emapper)

  output:
  tuple val(sample_id), path("eggmapper", type:"dir")
  
  script:
  """
  python ${projectDir}/bin/EggnogParser.py -n ${dir_emapper}/${sample_id}.emapper.annotations -f ${dir_emapper}/predict.faa -o eggmapper
  Rscript ${projectDir}/bin/draw_COG.R eggmapper/COG/all.COG.class.xls eggmapper/COG
  python ${projectDir}/bin/GO_anno.py -i eggmapper/${sample_id}.emapper.annotations -db ${params.db.GO_obo}
  Rscript ${projectDir}/bin/draw_GO.R eggmapper/GO/GO_anno_stats.xls eggmapper/GO
  python ${projectDir}/bin/KEGG_anno.py -i eggmapper/${sample_id}.emapper.annotations -db ${params.db.KEGG_dir}
  Rscript ${projectDir}/bin/draw_KEGG.R eggmapper/KEGG/KEGG_anno.stat.txt eggmapper/KEGG
  """
}
