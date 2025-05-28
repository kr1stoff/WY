#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.outdir = "./WholeGenomeAnalysis"
params.help = false
params.biwf = "batassy"
params.threads = 16
// params.csv


// 命令行参数
if (params.help) {
  log.info """
ONT+ILMN 全基因组无参组装分析流程
=============================================
用法:
  nextflow run ONTwholeGenome/nf-wholegenome
参数:
  * --csv     输入Samplesheet文件, .csv格式. ['Sample_ID','ILMN_R1','ILMN_R2','ONT_Read']
  * --outdir  输出目录. Default [${params.outdir}]
  * --biwf    生信工作流(Bioinfomatics Workflow)名称.
  * --help    打印该帮助文档.

生信工作流 (默认: batassy):
  batassy ::  细菌组装流程
  virassy ::  病毒组装流程
  """
  exit 0
}


// 导入工作流
include { BATASSY } from './workflows/batassy.nf'
include { VIRASSY } from './workflows/virassy.nf'


// channel: [sample_id, read1, read2, readlone]
Channel
    .fromPath( "${params.csv}" )
    .splitCsv()
    .map { row -> tuple(row[0], file(row[1]), file(row[2]), file(row[3])) }
    .set { ch_fastq }


// 主流程
workflow {
  if ( params.biwf == "virassy") {
    // 病毒组装流程
    VIRASSY(ch_fastq)
  }
  else {
    // 细菌组装流程
    BATASSY(ch_fastq)
  }
}