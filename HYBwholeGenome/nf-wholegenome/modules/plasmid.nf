process PLASMID_BLAST {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/plasmid", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(dir_predict)

  output:
  tuple val(sample_id), path("plasmid.blast")

  script:
  """
  blastn -num_threads ${params.threads} -evalue 1e-5 -max_target_seqs 5000 \
    -query ${dir_predict}/predict.ffn \
    -db $params.db.PLASMID_blastdb \
    -perc_identity 90 \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus' \
    -out plasmid.blast
  """
}


process PLASMID_FILTER {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/plasmid", mode: "copy"

  input:
  tuple val(sample_id), path(plasmid_blastout)

  output:
  tuple val(sample_id), path("plasmid.tsv")

  script:
  """
  #!/usr/bin/env python

  import pandas as pd

  def parse_tab(file):
      df = pd.read_csv(file,sep='\\t',encoding='utf-8',header=0)
      df.fillna('-',inplace=True)
      return df

  def parse_blast(file):
      df = pd.read_csv(file,sep='\\t',encoding='utf-8',header=None)
      columns = ['qid','sid','identity','length' ,'mismatch' ,'gapopen' ,
                'qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore' ,
                'qlen' ,'slen' ,'qcovs' ,'qcovhsp' ,'qcovus']
      df.columns = columns
      filter_df = df.groupby('sid').head(1)  # 比对到同一个参考只保留1条记录
      filter_df = filter_df.drop_duplicates(subset=['qid','qstart', 'qend'], keep='first')
      return filter_df

  tab_df = parse_tab("$params.db.PLASMID_annotbl")
  blast_df = parse_blast("$plasmid_blastout")

  df = pd.merge(blast_df,tab_df,how='left',left_on='sid',right_on='NUCCORE_ACC')
  df = df[df['NUCCORE_Completeness']=='complete']
  
  df = df.drop(columns=['NUCCORE_ACC','NUCCORE_Completeness'])
  out = "plasmid.tsv"
  df.to_csv(out,sep='\\t',encoding='utf-8',index=False)
  """
}