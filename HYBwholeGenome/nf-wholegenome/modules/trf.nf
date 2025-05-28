process TRF {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred/repeat", mode: "copy"
  conda params.conda_env.denovo

  input:
  tuple val(sample_id), path(fasta_asm)
  
  output:
  tuple val(sample_id), path("TandemRepeat.dat")
  
  script:
  """
  trf $fasta_asm 2 7 7 80 10 50 500 -h -m -ngs > TandemRepeat.dat
  """
}


process TRF_PARSE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.pred/repeat", mode: "copy"

  input:
  tuple val(sample_id), path(trf_dat)

  output:
  tuple val(sample_id), path("TandemRepeat.tsv")

  script:
  """
  #!/usr/bin/env python
  import re
  import pandas as pd
  
  def extract_number(s):
      match = re.search(r'\\d+', s)
      return int(match.group()) if match else 0

  def parse_Trf_Dat(file, out):
      df = pd.DataFrame()
      with open(file,'r',encoding='utf-8') as f:
          for line in f:
              if not line or line.startswith('#'): continue
              if line.startswith('@'):
                  name = line.strip().replace("@", "")
              else:
                  # TRF Start End Period Size Copy Number Consensus Size Percent 
                  # Matches Percent Indels Score A C G T Entropy(0-2) left reight
                  res = re.split(r'\\s',line)[0:13]
                  res = [name] + res
                  df = pd.concat([df,pd.DataFrame(res).T],axis=0,join='outer')
      columns = ['ID','Start','End','Period Size','Copy Number','Consensus Size',
              'Percent Matches','Percent Indels','Score','A','G','C','T','Entropy(0-2)']
      df.columns = columns
      df.sort_values('ID', key=lambda x: x.apply(extract_number),inplace=True)
      df.to_csv(out,sep='\\t',encoding='utf-8',index=False)

  parse_Trf_Dat("$trf_dat", "TandemRepeat.tsv")
  """
}