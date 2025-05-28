process ASSEMBLY_DEPTH {
  tag "$sample_id"
  publishDir "${params.outdir}/${sample_id}/2.asm/depth", mode: "copy"
  conda params.conda_env.basic

  input:
  tuple val(sample_id), path(fasta_asm)
  tuple val(sample_id), path(fastq)

  output:
  tuple val(sample_id), path("sort.bam.depth"), emit: depth_asm
  tuple val(sample_id), path("sort.bam"), emit: bam_asm
  
  script:
  """
  minimap2 -d ${fasta_asm}.mmi $fasta_asm 
  samtools faidx $fasta_asm
  minimap2 -t ${params.threads} -a ${fasta_asm}.mmi $fastq | samtools view -@ 16 -bS | samtools sort -@ 16 > sort.bam
  samtools depth -a sort.bam > sort.bam.depth
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
