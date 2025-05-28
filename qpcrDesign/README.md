# qpcrDesign
巢式PCR设计工具箱. 包含qPCR单个和批量设计,预扩增引物设计,引物简并转换


## 设计思路
- qpcr  
  基于`primer3`的qPCR设计脚本, 替代Windows版`PrimerExpress3.0.1`实现自动化引物&探针设计
  1. 输入FASTA都格式化成大写字母  
  2. 内引物设计, 使用`primer3`
  3. 内引物解析, 自编脚本

- qpcr-batch  
  批量版qpcr
  1. 从BED获取FASTA, 并拆分
  2. 使用`qpcr`批量跑, 初次设计使用70-110扩增产物长度阈值
  3. 使用`qpcr`批量跑, 针对70-110阈值不成功的区域再次设计,使用50-150扩增产物长度阈值

- prepcr  
  预扩增和tNGS设计方法一致, 使用`NGS-PrimerPlex`设计

- translate-degenerate-base
  简并碱基转化
  1. 输入巢式PCR经过引物评估后的最终表格  
  2. 输入本地数据库单物种基因组序列库 `all.fna`
  3. 提取内外引物探针序列
  4. `seqkit amplicon` 提取 `all.fna` 扩增子序列
  5. 检查扩增子中的简并和低频SNP  
    5.1. 获取本地单物种基因组条数和参考基因组长度, 分开 `ref` 和 `others` 序列  
    5.2. `mafft` 多序列比对, `--6merpair --keeplength --addfragments` 保持所有序列等长  
    5.3. 使用 字典中套 `ATGC` 数组, 记录每个碱基每个位置的频数   
    5.4. 根据简并多态性阈值来标记简并/小写/原始三种, 默认 `snp_ratio>5%`   
    5.5. 输出 `REF/A/T/G/C` 屏蔽保守的矩阵文件  
    5.6. 输出 `ref` 和 `degenerate` 的序列文件  
    5.7. 


## 参数调整
- 230529 内引物设计Primer3
  ```
  PRIMER_MAX_END_GC=3
  PRIMER_MAX_POLY_X=4
  PRIMER_INTERNAL_MAX_POLY_X=4
  PRIMER_GC_CLAMP=1
  ```
- 231115 内引物设计Primer3
  ```
  PRIMER_OPT_TM=64
  PRIMER_MIN_TM=62
  PRIMER_MAX_TM=67
  PRIMER_INTERNAL_OPT_TM=74
  PRIMER_INTERNAL_MIN_TM=71
  PRIMER_INTERNAL_MAX_TM=76

  PRIMER_SALT_CORRECTIONS=1
  PRIMER_TM_FORMULA=1
  PRIMER_SALT_MONOVALENT=150.0
  PRIMER_SALT_DIVALENT=1.5
  PRIMER_DNTP_CONC=0.6
  PRIMER_DNA_CONC=50.0
  PRIMER_INTERNAL_SALT_MONOVALENT=150.0
  PRIMER_INTERNAL_SALT_DIVALENT=1.5
  PRIMER_INTERNAL_DNTP_CONC=0.6
  PRIMER_INTERNAL_DNA_CONC=50.0
  ```


## 更新
- [x] 230915 引物评估转到这里  
- [x] 230710 巢式PCR外引物设计, 使用`NGS-PrimerPlex`, 失败可重试若干次替代手动检查log
- [x] 230704 `primer3_core`参数修复, OLD VERSION -> NEW VERSION。 https://primer3.org/manual.html
- [x] 230529 `primer3_core`部分更新, 使用Primer3Plus固定参数提高设计敏感性, 仅`PRODUCT_SIZE_RANGE`可变
- [x] 230519 批量脚本, Name列和BED关联
