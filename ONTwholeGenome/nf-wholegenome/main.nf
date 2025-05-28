#!/usr/bin/env nextflow

nextflow.enable.dsl=2
params.help = false

// 命令行参数
if (params.help) {
  log.info """
Nanopore 全基因组组装分析流程
=============================================
用法:
  nextflow run ONTwholeGenome/nf-wholegenome
参数:
  * --fastq   输入FASTQ文件, 样本名从该参数提取. e.g. ./.rawdata/*.fastq
  * --outdir  输出目录. Default [${params.outdir}]
  * --biwf    生信工作流(Bioinfomatics Workflow)名称.
  * --ref     参考基因组路径.
  * --help    打印该帮助文档.

生信工作流 (默认: batassy):
  batassy ::  细菌组装流程
  virassy ::  病毒组装流程
  virmap  ::  病毒有参拼接流程
  """
  exit 0
}

// parameters
params.outdir = "./WholeGenomeAnalysis"
params.biwf = "batassy"
params.threads = 16
// params.fastq
params.ref = '.'
def fileRef = new File(params.ref)
params.fileRefAbso = fileRef.absolutePath


// 导入工作流
include { BATASSY } from './workflows/batassy.nf'
include { VIRASSY } from './workflows/virassy.nf'
include { VIRMAP } from './workflows/virmap.nf'


// 原始数据
Channel.fromFilePairs(params.fastq, size: 1, checkIfExists: true).set{ ch_fastq }


// 每个subworkflow都写一个COMPLETE文件, 分析流程都跑完再生成报告
// 主流程
workflow {
  if ( params.biwf == "virassy") {
    // 病毒组装流程
    VIRASSY(ch_fastq)
  }
  else if ( params.biwf == "virmap") {
    // 病毒有参拼接流程
    assert fileRef.exists()
    assert fileRef.name != '.'
    VIRMAP(ch_fastq)
  }
  else {
    // 细菌组装流程
    BATASSY(ch_fastq)
  }
}