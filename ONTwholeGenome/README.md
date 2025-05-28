# ONTwholeGenome

三代细菌/病毒全基因组分析流程

## 使用

- 命令行参数

  ```
  usage: main.py [-h] -i INYAML [-o ANALYSIS] [-w {batassy,virassy,virmap}]

  ONT平台全基因组测序分析流程.

  optional arguments:
  -h, --help        show this help message and exit
  -i INYAML, --inyaml INYAML
                  输入YAML文件, SOLAR生成.
  -o ANALYSIS, --analysis ANALYSIS
                  结果分析目录, 在该目录下层新建文件夹保存结果, 比如 'analysis/task_id'. (default: .)
  -w {batassy,virassy,virmap}, --biwf {batassy,virassy,virmap}
                  生信工作流(Bioinfomatics Workflow)名称. batassy:细菌组装,virassy:病毒组装,virmap:病毒有参. (default: batassy)
  ```

- 实例

  ```bash
  # 细菌组装
  python ONTwholeGenome/main.py -i ONTwholeGenome/template/bacassy.yaml -o result

  # 病毒组装
  python ONTwholeGenome/main.py -w virassy -i ONTwholeGenome/template/virassy.yaml -o result

  # 病毒有参拼接
  python ONTwholeGenome/main.py -w virmap -i ONTwholeGenome/template/virmap.yaml -o result
  ```

## 测试数据

- 细菌 https://journals.asm.org/doi/10.1128/mra.01129-22
- 病毒 https://www.sciencedirect.com/science/article/pii/S0166093421001191

## 运行步骤

软件运行流程

1. 解析 SOLAR MySQL YAML 文件
2. 创建部分目录
3. 链接原始数据 (FASTQ 或 zcat FASTQ_DIR/\*)
4. 运行 Nextflow 主脚本
5. 拷贝结果到 Upload 目录
6. 生成报告
7. 压缩结果目录

## Nextflow

### 细菌无参

1. 质控
   - NanoStat
   - NanoPlot
   - NanoFilt
2. 组装
   - Jellyfish
   - Genomescope2
   - FLYE
   - QUAST
   - CheckM
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
   - CANU
   - MEDAKA_CONSENSUS
   - CHECKV
   - QUAST
3. 基因预测
   - Prokka
4. 基因注释
   - Swiss-Prot
   - KEGG

### 病毒有参

1. 质控 (同上)
2. 比对
   - Minimap2
3. 变异
   - LoFreq
4. 拼接
   - Bedtools

## 更新

- [x] 240523 内容升级
  - 组装信息统计表加成环状态信息
- [x] 240311 线程数 BUG 修复
  - nextflow.config 中 cpus 设置超过当前服务器最大线程数报错.
- [x] 240227 流程更新
  - SOLAR 前端联调
  - virmap BUG 修复. 由于 `ch_ref` 在比对部分只会运行一个样本, 删除 `ch_ref` 改成 `params.ref` 传递参数
- [x] 240105 内容更新
  - 串联重复预测
  - 毒力注释优化
  - 耐药部分新增 ResFinder
  - 新增质粒注释
