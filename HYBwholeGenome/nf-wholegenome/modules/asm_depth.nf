process ASSEMBLY_DEPTH {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/depth", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(fasta_asm)
  tuple val(sample_id), path(fastq)
  tuple val(sample_id), path(read1), path(read2)

  output:
  tuple val(sample_id), path("merged.bam.depth"), emit: depth_asm
  tuple val(sample_id), path("merged.bam"), emit: bam_asm
  
  script:
  """
  # 索引
  minimap2 -d ${fasta_asm}.mmi $fasta_asm
  bwa index -a bwtsw $fasta_asm
  samtools faidx $fasta_asm
  # 比对
  minimap2 -t ${params.threads} -a ${fasta_asm}.mmi $fastq | samtools view -@ ${params.threads} -bS \
    | samtools sort -@ ${params.threads} > sort.long.bam
  bwa mem -t ${params.threads} $fasta_asm $read1 $read2 | samtools view -@ ${params.threads} -bS \
    | samtools sort -@ ${params.threads} > sort.short.bam
  # 合并BAM
  printf '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:Nanopore\\n@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA\\n' > rg.txt
  samtools merge -rh rg.txt merged.bam sort.long.bam sort.short.bam
  # 深度
  samtools depth -a merged.bam > merged.bam.depth
  """
}


process GC_DEPTH_CONSISTENCY {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/depth", mode: "copy"

  input:
  tuple val(sample_id), path(file_depth)
  tuple val(sample_id), path(fasta_asm)

  output:
  tuple val(sample_id), 
        path("consistency.txt"), 
        path("depth_base.stat.depth_GC.png"), 
        path("consistency.all.txt")
  
  script:
  """
  python ${projectDir}/bin/depth_base_stat.py -l 2000 -g $fasta_asm -d $file_depth -s depth_base.stat
  Rscript ${projectDir}/bin/depth_GC_plot.r -i depth_base.stat -o depth_base.stat.depth_GC
  convert depth_base.stat.depth_GC.pdf depth_base.stat.depth_GC.png
  python ${projectDir}/bin/count_depth.py -i $file_depth -ref $fasta_asm -o consistency.txt -a consistency.all.txt
  """
}
