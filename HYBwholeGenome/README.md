# HYBwholeGenome

ONT+ILMN 混装流程

## 命令行

- 命令行参数

  ```
  usage: main.py [-h] -i INYAML [-o ANALYSIS] [-w {batassy,virassy}]

  ONT+ILMN混装全基因组测序分析流程.

  optional arguments:
   -h, --help         show this help message and exit
   -i INYAML, --inyaml INYAML
                    输入YAML文件, SOLAR生成.
   -o ANALYSIS, --analysis ANALYSIS
                    结果分析目录, 在该目录下层新建文件夹保存结果, 比如 'analysis/task_id'. (default: .)
   -w {batassy,virassy}, --biwf {batassy,virassy}
                    生信工作流(Bioinfomatics Workflow)名称. batassy:细菌组装,virassy:病毒组装. (default: batassy)=
  ```

- 实例

  ```bash
  # 细菌无参
  python HYBwholeGenome/main.py -i HYBwholeGenome/template/bacassy.yaml -o result

  # 病毒无参
  python HYBwholeGenome/main.py -w virassy -i HYBwholeGenome/template/virassy.yaml -o result
  ```

## 测试数据

- 细菌 https://journals.asm.org/doi/10.1128/mra.01129-22
- 病毒 https://www.mdpi.com/1999-4915/15/6/1248

## 运行步骤

软件运行流程

1. 解析 SOLAR MySQL YAML 文件
2. 创建部分目录
3. 链接原始数据 (FASTQ 或 zcat FASTQ_DIR/\*)
4. 运行 Nextflow 主脚本
5. 拷贝结果到 Upload 目录
6. 生成报告
7. 压缩结果目录

## nextflow

### 细菌

1. 数据质控
   - NanoStat
   - NanoPlot
   - NanoFilt
   - FastQC
   - fastp
2. 组装
   - Unicycler
   - QUAST
   - CHECKM
3. 基因预测
   - Prokka
   - Phispy
   - Island
4. 基因注释
   - VFDB
   - CARD
   - EggNOG-mapper
   - Cazy
   - Swiss-Prot

### 病毒无参

1. 质控 (同上)
2. 组装
   - unicycler
   - QUAST
   - CHECKM
3. 基因预测
   - Prokka
4. 基因注释
   - Swiss-Prot
   - KEGG

## 更新

- [x] 240523 内容升级
  - 组装信息统计表加成环状态信息
- [x] 240311 线程数 BUG 修复
  - nextflow.config 中 cpus 设置超过当前服务器最大线程数报错.
- [x] 240301 流程更新
  - SOLAR 前端联调
