# tngsAnalysis
WY扩增子平台结果分析部分.  

## 用法
```
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  WY扩增子平台结果统计.

Options:
  -h, --help  Show this message and exit.

Commands:
  dimer-assessment  引物二聚体评估.
```

## 设计思路
### 引物二聚体评估 dimer-assessment
- 单端测序
  1. `seqtk seq -A` 输入FASTQ转FASTA;  
  2. 引物序列建库，`blast`比对测序序列到引物上;  
  3. 整理字典1, 以qseqid为主key, 将所有qseqid比对到的位置写入该qseqid下
  4. 整理字典2, 根据字典1将qseqid下面的所有sseqid以及位置信息连接起来，整理扩增子比对位置的出现位置的统计, 一条扩增子只比对到一个引物默认不是二聚体
  5. 绘图二聚体组用饼图和柱形图统计，最高的几个二聚体组柱形图统计最高出现的引物位置
- 双端测序  
  1. 在单端测序的基础上进行部分修改  
  2. 使用`Pandaseq`合并双端FASTQ，成对/不成对的都要输出，然后合并成1个FASTA  
  3. 主题思路和单端比对统计一样。统计部分如果一个扩增子只有两条引物，且这两条引物是一对，则认为该条扩增子不属于二聚体
- 结果  
  - primer_groups_count.tsv  
    ```
    group                       count
    159_R-295_R                 1912
    251_L-243_L                 1742
    283_R-203_R-203_R-215_L     1733
    ```
  - primer_groups_position_count.tsv
    ```
    group                    position                                           count
    307_R-307_R              89-113,1-25;237-261,1-25                           847
    327_L-327_L-191_R-191_R  37-58,1-22;186-207,1-22;47-56,10-19;196-205,10-19  714
    327_L-191_L-191_L        284-302,1-19;236-245,24-15;107-116,24-15           691
    ```
