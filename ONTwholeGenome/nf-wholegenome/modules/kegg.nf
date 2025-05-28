process KEGG {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/KEGG", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(dir_predict)
  
  output:
  tuple val(sample_id), path("detail.tsv")
  
  script:
  """
  exec_annotation -f detail-tsv -E 1e-3 --profile ${params.db.KOfam_dir}/profiles \
    --ko-list ${params.db.KOfam_dir}/ko_list --cpu ${params.threads} -o detail.tsv ${dir_predict}/predict.faa
  """
}

process KEGG_ANNO_STATS {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/KEGG", mode: "copy"

  input:
  tuple val(sample_id), path(kegg_detail)

  output:
  tuple val(sample_id), path("KEGG_anno.stat.tsv"), path("KEGG_anno_stats.png")

  script:
  """
  #!/usr/bin/env python
  import pandas as pd
  import seaborn as sns
  import matplotlib.pyplot as plt

  # KO LevelA LevelB 对照表
  kegg_kab = "${params.db.KOfam_dir}/KO_LevelA_LevelB.tsv"

  # 统计表
  df_kd = pd.read_table("$kegg_detail", delimiter='\\t', usecols=['gene name', 'KO', 'KO definition'])
  df_kd.drop(index=0, inplace=True)
  df_kab = pd.read_table(kegg_kab, sep='\\t')

  df_kegg_anno =  pd.merge(df_kd, df_kab, how='left', on='KO')
  df_kegg_anno.to_csv('KEGG_anno.stat.tsv', index=False, sep='\\t')


  # level A-B 统计图
  df_count_lvl_ab = df_kegg_anno[['gene name', 'LevelA', 'LevelB']].groupby(['LevelA', 'LevelB']).count().reset_index()
  # 创建图形
  sns.set(style="whitegrid")
  fig, ax = plt.subplots(figsize=(12, 7))
  # 绘制水平条形图
  barplot = sns.barplot(data=df_count_lvl_ab, y='LevelB', x='gene name', hue='LevelA', orient='h')
  # 移除图例的标题（标签）
  ax.legend_.remove()
  # 设置图例的位置（在图形之外）
  plt.legend(loc='best', bbox_to_anchor=(1, 0.75), title=None)
  # 添加坐标轴标签
  plt.xlabel('Gene Count')
  plt.ylabel('')
  # 显示图形
  plt.show()
  # 保存图片
  fig.savefig('KEGG_anno_stats.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
  """
}