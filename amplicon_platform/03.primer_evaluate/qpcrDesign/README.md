# qpcrDesign
- 基于`primer3`的qPCR设计脚本, 替代Windows版`PrimerExpress3.0.1`实现自动化引物&探针设计.
- 预扩增和tNGS设计方法一致, 使用`NGS-PrimerPlex`设计

## 设计思路
- qpcr  
  1. 输入FASTA都格式化成大写字母  
  2. 内引物设计, 使用`primer3`
  3. 内引物解析, 自编脚本

- qpcr-batch  
  1. 从BED获取FASTA, 并拆分.
  2. 使用`qpcr`批量跑, 初次设计使用70-110扩增产物长度阈值.
  3. 使用`qpcr`批量跑, 针对70-110阈值不成功的区域再次设计,使用50-150扩增产物长度阈值.

## 参数调整
- 230529参数调整  
  ```
  PRIMER_MAX_END_GC=3(5)
  PRIMER_MAX_POLY_X=4(5)
  PRIMER_INTERNAL_MAX_POLY_X=4(5)
  PRIMER_GC_CLAMP=1(0)
  ```

## 更新
- [x] 230710 巢式PCR外引物设计, 使用`NGS-PrimerPlex`, 失败可重试若干次替代手动检查log
- [x] 230704 `primer3_core`参数修复, OLD VERSION -> NEW VERSION. https://primer3.org/manual.html
- [x] 230529 `primer3_core`部分更新, 使用Primer3Plus固定参数提高设计敏感性, 仅`PRODUCT_SIZE_RANGE`可变
- [x] 230519 批量脚本, Name列和BED关联

## TODO
