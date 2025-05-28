process LOFREQ {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.muts", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(bam)
  tuple val(sample_id), path(ref)

  output:
  tuple val(sample_id), path("muts.vcf")

  script:
  """
  lofreq call --force-overwrite -q 20 -Q 20 -C 20 -f $ref -o muts.vcf $bam
  """
}


process EXTRECT_LOFREQ {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/3.muts", mode: "copy"

  input:
  tuple val(sample_id), path(vcf)

  output:
  tuple val(sample_id), path("muts.tsv")
  
  script:
  """
  #!/usr/bin/env python
  import re

  with open("$vcf") as f, open("muts.tsv", "wt") as g:
      g.write("CHROM\\tPOS\\tREF\\tALT\\tDEPTH\\tFREQ\\n")
      for line in f:
          if line.startswith('#'):
              continue
          lns = line.strip().split('\\t')
          chrom, pos, _, ref, alt = lns[:5]
          info = lns[7]
          res = re.findall("^DP=(\\d+);AF=(.*?);", info)
          if res:
              dp = res[0][0]
              af = "{:.2%}".format(float(res[0][1]))
              g.write("\\t".join([chrom, pos, ref, alt, dp, af]) + "\\n")
  """
}
