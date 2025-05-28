#!/usr/bin/env perl

use lib "/sdbb/share/lib/GDHR";
use strict;
use GDHR;
use Utils;
use File::Basename;

# 命令函参数
my $prog = basename($0);
my $usage = "Usage:\n\t$prog <upldsmp> <sample_name>\n";
die $usage if $#ARGV != 1;
my $upldsmp = $ARGV[0];
my $sample = $ARGV[1];

# 标题
my $report;
$report = GDHR->new(-outdir  => $upldsmp,
                    -pipe    => "新冠污水监测分析报告",
                    -nonlazy => 1);
# 项目概述
my $section;
$section = $report->section(id => "introduction");
$section->menu("项目概述", -icon => "source/icons/gaishu.png");
$section->desc("有研究表明, 相同地区污水检出新冠重点关注的分型时间会早于临床样本1-2周, 因此对污水分析可用于人群中新冠流行的监测。随着时间的推移监测人群水平的新冠变异和进化, 对新冠疫苗, 疗法和诊断的有效性至关重要。");


# 数据质控
my $section;
$section = $report->section(id => "dataQC");
$section->menu("数据质控", -icon => "source/icons/cexuzhikong.png");
FileLink($section, "原始数据质控结果请点击 <link=./1.qc>");

## 过滤前fastqc
$section->submenu("测序数据质量分析", -icon => "source/icons/guolv.png");
$section->desc("过滤前碱基质量值图, 展示原始下机序列上每个位置碱基的平均质量值, 结果如下图所示: ");
if(-e "$upldsmp/1.qc/$sample.2_before_per_base_quality.png"){ # 双端测序两个图
    my @name = [ "imgs2html2" ];
    $section->img2html2(
        -file1 => "1.qc/$sample.1_before_per_base_quality.png",
        -name1 => "双端测序序列1碱基质量图",
        -file2 => "1.qc/$sample.2_before_per_base_quality.png",
        -name2 => "双端测序序列2碱基质量图",
        -names => "img2html2");
}else{ # 单端测序一个图
    $section->img2html(
        -file => "1.qc/$sample\_before_per_base_quality.png",
        -name => "单端测序序列碱基质量图");
};
## 过滤后fastqc
$section->desc("过滤后碱基质量值图, 展示过滤后序列上每个位置碱基的平均质量值, 结果如下图所示: ");
if(-e "$upldsmp/1.qc/$sample.2_after_per_base_quality.png"){ # 双端测序两个图
    my @name = [ "imgs2html2" ];
    $section->img2html2(
        -file1 => "1.qc/$sample.1_after_per_base_quality.png",
        -name1 => "双端测序序列1碱基质量图",
        -file2 => "1.qc/$sample.2_after_per_base_quality.png",
        -name2 => "双端测序序列2碱基质量图",
        -names => "img2html2");
}else{ # 单端测序一个图
    $section->img2html(
        -file => "1.qc/$sample\_after_per_base_quality.png",
        -name => "单端测序序列碱基质量图");
};

## 数据过滤前后的对比结果
$section->submenu("数据统计", -icon => "source/icons/guolv.png");
$section->tsv2html( # 添加一个表格
    -file   => "$upldsmp/1.qc/$sample.basic.stat.txt",
    -top    => 100,
    -name   => "过滤指标统计",
    -header => 1);


# 结果分析
my $section;
$section = $report->section(id => "results");
$section->menu("结果分析", -icon => "source/icons/fenxi.png");

## 比对结果
$section->submenu("比对结果", -icon => "source/icons/bidui.png");
FileLink($section, "比对结果请点击 <link=./2.align>");
$section->desc("将过滤后 FASTQ 比对回参考基因组, 对比对结果进行统计。结果如下表所示: ");
$section->tsv2html(
    -file   => "$upldsmp/6.demix/$sample.coverage.txt",
    -top    => 100,
    -name   => "基因比对覆盖度与深度表",
    -header => 1);
$section->desc("结果说明<br>平均深度: 基因组平均覆盖深度, 覆盖度: 基因组覆盖度, 深度≥10x: 覆盖深度达到10层的区域在基因组中的占比, 均一性: 基因组覆盖深度的一致性指数。");
$section->desc("全基因组覆盖图如下图所示: ");
$section->img2html(
    -file => "6.demix/$sample.coverage.png",
    -name => "全基因组覆盖度图");
$section->desc("结果说明<br>横坐标为参考基因组位置区间, 纵坐标为该区间的平均测序覆盖深度。");

## 分型丰度
$section->submenu("分型丰度", -icon => "source/icons/xinguan.png");
$section->desc("污水样本中各分型丰度, 结果如下表所示: ");
$section->tsv2html(
    -file   => "$upldsmp/6.demix/$sample.abundance.txt",
    -top    => 100,
    -name   => "分型丰度结果表",
    -header => 1);
if(-e "$upldsmp/6.demix/$sample.abundance_user_predict.txt"){
    $section->desc("预测的分型, 结果如下表所示: ");
    $section->tsv2html(
        -file   => "$upldsmp/6.demix/$sample.abundance_user_predict.txt",
        -top    => 100,
        -name   => "预测分型丰度结果表",
        -header => 1);
}
$section->desc("污水样本中各分型丰度如下图所示: ");
$section->img2html(
    -file => "6.demix/$sample.pie.png",
    -name => "全基因组覆盖度图");

## 变异结果
$section->submenu("变异信息", -icon => "source/icons/bianyi.png");
$section->desc("污水样本中的位点变异信息及对应分型, 结果如下表所示: ");
$section->tsv2html(
    -file   => "$upldsmp/3.muts/$sample.display.tsv",
    -top    => 100,
    -name   => "变异结果表",
    -header => 1);
FileLink($section, "变异信息结果请点击 <link=3.muts/$sample.display.tsv>");

# 生成HTML报告
$report->write();
