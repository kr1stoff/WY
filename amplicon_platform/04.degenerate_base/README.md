# tngsPangenome
WY扩增子平台引物设计部分泛基因组分析.  


## 用法
```
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  WY扩增子平台引物设计泛基因组.

Options:
  -h, --help  Show this message and exit.

Commands:
  get-accession-list         ...
  get-conserved-from-roary   使用roary结果和gff文件坐标, 汇总表格.
  get-coregene-fa            获取同物种下各基因组特定核心基因序列, 写入一个FASTA文件中.
  get-extra-roary-stats      roary额外的统计
  get-predict-schedule       获取预测进度表.
  run-batch-prokka           批量跑prokka, 输入为物种目录包含多个基因组FASTA, 输出为GFF物种目录.
  translate-degenerate-base  简并碱基转换.
```


## 更新
- [ ] 230809 简并碱基转化软件
- [x] 221206 批量跑prokka脚本  
- 230104 BUG修复  
  保留参考基因组的收录号. 基因组超过两千会筛选收录号, 存在参考基因组收录号被过滤掉的情况.
- 230103 流程升级  
  - 小程序收录号列表文件, 软连接这些收录号的fna  
  - 小程序预测进度表  
  - 子程序名优化  
