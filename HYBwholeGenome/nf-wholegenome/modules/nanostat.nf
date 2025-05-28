process NANOSTAT_PRE {
  // process 标记
  tag "$sample_id"
  // 输出目录
  publishDir "${params.outdir}/${sample_id}/1.qc/nanostat", mode: "copy"
  // 使用现有的CONDA环境
  conda params.conda_env.nano
  
  input:
  tuple val(sample_id), path(read1), path(read2), path(fastq)

  output: 
  tuple val(sample_id), path("${sample_id}.raw.nanostat")

  script:
  """
  NanoStat --fastq $fastq > ${sample_id}.raw.nanostat
  """
}


process NANOSTAT_POST {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc/nanostat", mode: "copy"
  conda params.conda_env.nano

  input:
  tuple val(sample_id), path(fastq)

  output: 
  tuple val(sample_id), path("${sample_id}.cln.nanostat")

  script:
  """
  NanoStat --fastq $fastq > ${sample_id}.cln.nanostat
  """
}

process EXTRACT_NANOSTAT {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/1.qc/nanostat", mode: "copy"

  input:
  tuple val(sample_id), path(nanostat_raw)
  tuple val(sample_id), path(nanostat_cln)

  output:
  tuple val(sample_id), path("${sample_id}.nanostat")

  script:
  """
  #!/usr/bin/env python
  from collections import OrderedDict

  def get_nanostat_dict(nanostat):
      with open(nanostat) as f:
          dicnnst = OrderedDict()
          for line in f:
              if ':' in line:
                  lst = line.strip().split(':')
                  dicnnst[lst[0]] = lst[1].strip()
      return dicnnst

  dicraw = get_nanostat_dict("${nanostat_raw}")
  diccln = get_nanostat_dict("${nanostat_cln}")
  with open("${sample_id}.nanostat", 'w') as g:
      g.write('过滤质控\\t过滤前\\t过滤后\\n')
      g.write(f"序列数\\t{dicraw['Number of reads']}\\t{diccln['Number of reads']}\\n")
      g.write(f"总碱基数\\t{dicraw['Total bases']}\\t{diccln['Total bases']}\\n")
      g.write(f"序列平均长度\\t{dicraw['Mean read length']}\\t{diccln['Mean read length']}\\n")
      g.write(f"序列长度N50\\t{dicraw['Read length N50']}\\t{diccln['Read length N50']}\\n")
      g.write(f"平均序列质量\\t{dicraw['Mean read quality']}\\t{diccln['Mean read quality']}\\n")
      g.write(f">Q7\\t{' '.join(dicraw['>Q7'].split(' ')[:2])}\\t{' '.join(diccln['>Q7'].split(' ')[:2])}\\n")
      g.write(f">Q12\\t{' '.join(dicraw['>Q12'].split(' ')[:2])}\\t{' '.join(diccln['>Q12'].split(' ')[:2])}\\n")
  """
}