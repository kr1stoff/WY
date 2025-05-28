# PathoGenomeIntegrator

病原基因组整合器，单个物种的多个基因组合并

## 思路

1. 各株系基因组生成假 fastq 数据
2. 比对 fastq 到参考基因组，获取 unmapped reads
3. 组装 unmapped reads ，合并参考序列生成 mixed libraries（整合基因组） 参考
4. 比对 fastq 到 mix libraries 参考

### 未指定参考基因组

如果软件参数未指定参考基因组登录号, 需要自行筛选.

1. checkM
    - ~~单拷贝标记物 (完整度和污染度, 结合了缺失,单拷贝,多拷贝计算的)~~
    - 完整度
    - 污染度
2. prokka
    - 计算假基因数量占比
      prokka.tsv 文件中, (gene rows count) - (CDS), 可以理解成假基因数量
3. seqtk comp
    - N 碱基占比
      结果中 "#4" 列是 N 碱基, 除以总长度就是 N 碱基比例
    - contig 数量
      根据 RefSeq assembly_summary_refseq 表格中, 所有统计到的基因组(大概 35 万+), 分类级别在"Chromosome"和"Complete Genome"级别 contig 数量大概在 4 以下(75%), "Contig"和"Scaffold"级别 contig 数量大概在 130 以上(75%)
