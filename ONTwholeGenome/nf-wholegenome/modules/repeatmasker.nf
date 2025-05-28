process REPEATMASKER {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred", mode: "copy"
  conda params.conda_env.denovo

  input:
  tuple val(sample_id), path(fasta_asm)
  
  output:
  tuple val(sample_id), path("repeat", type:"dir")

  when:
  params.biwf == "batassy"

  script:
  """
  RepeatMasker -e ncbi $fasta_asm -dir repeat -pa ${params.threads} -s -nolow -gff -lib ${params.db.RepeatMasker_lib_dir}/RepeatMasker.lib
  """
}


process REPEATMASKER_PARSE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred/repeat", mode: "copy"

  input:
  tuple val(sample_id), path(dir_repeat)

  output:
  tuple val(sample_id), path("INEs.tsv")

  when:
  params.biwf == "batassy"

  script:
  """
  #!/usr/bin/env python
  import re
  import pandas as pd

  def extract_number(s):
      match = re.search(r'\\d+', s)
      return int(match.group()) if match else 0

  def parse_file(file):
      res = []
      with open(file,'r',encoding='utf-8') as f:
          for line in f:
              if not line or line.startswith('#'): continue
              args = line.strip().split("\\t")
              queryid,start,end,strand,kind = args[0],args[3],args[4],args[6],args[8]
              queryid = queryid.replace("@", "")
              kind = re.findall(r":(.*?)\\"" ,kind)[0]
              res.append([queryid,start,end,strand,kind])
      df = pd.DataFrame(res,columns=['qid','qstart','qend','strand','class'])
      df.sort_values('qid', key=lambda x: x.apply(extract_number),inplace=True)
      return df
      
  def parse_anno(anno):
      df = pd.read_csv(anno,sep='\\t',encoding='utf-8')
      return df

  df = parse_file("${dir_repeat}/assembly.fasta.out.gff")
  anno_df = parse_anno("${params.db.RepeatMasker_lib_dir}/RepeatMasker.anno.tab")
  merge_df = pd.merge(df,anno_df,left_on='class',right_on='ID',how='inner')
  merge_df = merge_df[['qid','qstart','qend','strand','class','Type','SubType','Species']]
  merge_df.to_csv("INEs.tsv",sep='\\t',encoding='utf-8',index=False) 
  """
}