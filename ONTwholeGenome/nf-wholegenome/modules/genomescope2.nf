process GENOMESCOPE2 {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm", mode: "copy"
  conda params.conda_env.asm

  input:
  tuple val(sample_id), path(jf_histo)

  output: 
  tuple val(sample_id), path("genomescope2/${sample_id}_summary.txt")

  script:
  """
  genomescope2 -i ${jf_histo} -o genomescope2 -n ${sample_id} -p 1 -k 21
  """
}


process GENOMESIZE1 {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/genomescope2", mode: "copy"

  input:
  tuple val(sample_id), path(genomescope2_sum)

  output:
  tuple val(sample_id), path("genome_size.txt")

  script:
  """
  #!/usr/bin/env python
  import re

  with open("${sample_id}_summary.txt") as f:
      for line in f:
          if "Genome Haploid Length" in line:
              lns = re.split("  +", line)
              min_len = lns[1].replace(" bp","").replace(",","")
              max_len = lns[2].replace(" bp","").replace(",","")    # "Inf bp" 报错
              if max_len.isdigit():
                  gnmlen = int(max_len)
              elif min_len.isdigit():
                  gnmlen = int(min_len)
              else:
                  gnmlen = 1000001
              if gnmlen >= 1000000:
                  lth_fmt = f"{(gnmlen/1000000):.1f}m"
              else:
                  lth_fmt = f"{(gnmlen/1000):.1f}k"
      with open("genome_size.txt", "w") as g:
          g.write(lth_fmt)
  """
}


process GENOMESIZE2 {
  tag "$sample_id"

  input:
  tuple val(sample_id), path(fl_genome_size)

  output:
  val(val_genome_size)

  exec:
  def file = new File("${params.outdir}/${sample_id}/2.asm/genomescope2/genome_size.txt")
  val_genome_size = file.text
}
