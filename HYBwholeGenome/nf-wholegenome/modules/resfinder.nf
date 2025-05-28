process RESFINDER {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno", mode: "copy"
  conda params.conda_env.resfinder

  input:
  tuple val(sample_id), path(dir_predict)

  output:
  tuple val(sample_id), path("ResFinder", type:"dir")

  script:
  """
  run_resfinder.py -ifa ${dir_predict}/predict.ffn --db_path_res ${params.db.Resfinder_dir} -s "Other" -o ResFinder \
    -acq -l 0.8 -t 0.8 --nanopore
  """
}


process RESFINDER_PARSE {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/4.anno", mode: "copy"

  input:
  tuple val(sample_id), path(dir_resfinder)
  
  output:
  tuple val(sample_id), path("ResFinder/Resfinder_anno.txt")
  
  script:
  """
  #!/usr/bin/env python

  import os
  import pandas as pd

  class ParseResfinder():
      def __init__(self,**kwargs):
          self.db_file = os.path.abspath(kwargs['db_file'])
          self.file = os.path.abspath(kwargs['infile'])
          self.flag = True if os.path.exists(self.file) else False
          self.outfile = os.path.abspath(kwargs['outfile'])

      def parse_file(self, file):
          df = pd.read_csv(file,sep='\\t',encoding='utf-8',comment='#')
          df.fillna('-',inplace=True)
          return df

      def merge_df(self, df1, df2, key):
          df = pd.merge(df1, df2, on=key, how='inner')
          return df

      def run(self):
          OUT = os.path.join(self.outfile)
          if self.flag:
              anno_df = self.parse_file(self.file)
              db_df = self.parse_file(self.db_file)
              db_df['Accession no.'] = db_df['Gene_accession no.'].str.extract(r'_.*?_(.*?)\$')
              merge_df = self.merge_df(anno_df,db_df,key='Accession no.')
              del merge_df['Notes']
              merge_df.to_csv(OUT,sep='\\t',encoding='utf-8',index=False)
          else:
              with open(OUT,'w',encoding='utf-8') as w:
                  w.write(f"Resistance gene\\tIdentity\\tAlignment Length/Gene Length\\tCoverage\\tPosition in reference\\tContig\\tPosition in contig\\tPhenotype_x\\tAccession no.\\tGene_accession no.\\tClass\\tPhenotype_y\\tPMID\\tMechanism of resistance\\tRequired_gene")

  project = ParseResfinder(
      infile = "${dir_resfinder}/ResFinder_results_tab.txt",
      db_file = "${params.db.Resfinder_dir}/phenotypes.txt",
      outfile = "ResFinder/Resfinder_anno.txt"
  )
  project.run()
  """
}