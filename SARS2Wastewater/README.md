# 新冠污水检测流程
## 安装  
创建snakemake环境  
软件依赖  
- snakemake
- pandas
- numpy
- matplotlib
- seaborn

本地需要安装的其他软件  
- perl
- fastqc
- fastp
- bowtie2
- ivar
- samtools
- bcftools
- pangolin
- freyja


## 数据库文件
- usher_barcodes.csv  
  https://github.com/andersen-lab/Freyja/blob/main/freyja/data/usher_barcodes.csv  


## 注意  
1. `freyja`创建新环境使用, 其中要注意`pandas=1.3.5`, 高版本的`pandas`已经移除`DataFrame.append`函数, **会报错!!**  


## UPDATE
- [x] [230423] 三代测序流程
- [x] [230417] 支持primer参数, 去引物bed文件输入  
- [x] [230417] 用户输入预测分型
- [x] [230417] 变异表格仅展示freyja优势物种
- [x] [230414] `freyja` 分型丰度表全展示  
