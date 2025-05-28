process VFDB_BLAST {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/VFDB", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(dir_predict)
  
  output:
  tuple val(sample_id), path("blast.out")

  script:
  """
  blastn -num_threads ${params.threads} -evalue 1e-10 -query ${dir_predict}/predict.ffn -db ${params.db.VFDB_blastdb} \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus' \
    -out blast.out
  """
}


process VFDB_PARSE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/VFDB", mode: "copy"

  input:
  tuple val(sample_id), path(blast_out)
  tuple val(sample_id), path(depth_predict)
  
  output:
  tuple val(sample_id), path("Virulence.tsv")
  
  script:
  """
  #!/usr/bin/env python

  import pandas as pd

  def parse_blast_out(file):
      # 解析blast结果
      columns = ['Contig','Gene ID','identity','length' ,'mismatch' ,'gapopen' ,
                'qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore' ,
                'qlen' ,'slen' ,'qcovs' ,'qcovhsp' ,'qcovus']
      df = pd.read_csv(file, sep='\\t', encoding='utf-8', header=None)
      df.columns = columns
      # df = df.groupby('Contig').head(1)  # 比对到同一个参考只保留1条记录
      df = df.drop_duplicates(subset=['Gene ID','qstart', 'qend'], keep='first')
      #条件筛选，覆盖度
      df = df.loc[(df["qcovs"] > 80)]
      return df

  def parse_tab(file):
      df = pd.read_csv(file, sep='\\t', encoding='utf-8', header=0)
      df.fillna('-',inplace=True)
      return df

  blast_df = parse_blast_out("$blast_out")
  ref_df = parse_tab("$params.db.VFDB_merge_tab")

  depcov_df = parse_tab("$depth_predict")
  depcov_df.columns = ["Contig", "Depth", "Coverage"]
  depcov_df["Depth"] = depcov_df['Depth'].apply(lambda x: '{:.2f}'.format(x))
  depcov_df["Coverage"] = depcov_df['Coverage'].apply(lambda x: '{:.2f}'.format(x))
  total_df = pd.merge(blast_df, depcov_df, how='left', on="Contig")
  total_df = pd.merge(total_df, ref_df, how='left', on="Gene ID")
  
  total_df = total_df.drop_duplicates(subset=['Related Genes','Virulence Factors'], keep='first')
  total_df.fillna('-',inplace=True)
  total_df.to_csv("Virulence.tsv", sep='\\t',encoding='utf-8',index=False)
  """
}


process VFDB_FASTA {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno/VFDB", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(virulence_gene)
  tuple val(sample_id), path(dir_predict)
  
  output:
  tuple val(sample_id), path("viru.fa"), path("viru_ref.fa")
  
  script:
  """
  awk 'NR>1' $virulence_gene | cut -f1 | sort | uniq > show_virulence_gene.id
  seqkit grep -f show_virulence_gene.id ${dir_predict}/predict.ffn > viru.fa
  awk 'NR>1' $virulence_gene | cut -f2 | sort | uniq > ref.id
  seqkit grep -f ref.id ${params.db.VFDB_setB_fa} > viru_ref.fa
  """
}