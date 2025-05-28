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
my $dir = $ARGV[0];
my $sample = $ARGV[1];

my $blank6 = blank(6);


# 标题
my $report;
$report = GDHR->new(-outdir => $dir,
    -pipe                   => "病毒全基因组组装分析报告",
    -nonlazy                => 1);


# 1 项目介绍
my $section;
$section = $report->section(id => "introduction");
$section->menu("项目概述");
$section->desc("$blank6 通过分离培养获得病毒株，进行建库测序，获得较为纯净的毒株序列。本项目对单毒株（样本编号为$sample）的Nanopore测序数据进行分析，通过输入参考基因组，进行比对拼接得到一致性序列。");


# 2 数据质控
my $section;
$section = $report->section(id => "stat_show");
$section->menu("数据质控",-icon => "./src/icon/2.raw_qc.png");
FileLink($section, "$blank6  全部原始数据质控结果请点击  <link=./1.qc>");

# 数据过滤介绍
$section->desc("$blank6 测序数据的产生是经过了 DNA 提取、建库、测序多个步骤的，这些步骤中产生的无效数据会对后续信息分析带来严重干扰，如建库长度的偏差，以及测序错误、低质量碱基、未测出的碱基（以 N 表示）等情况，须通过一些手段将上述无效数据过滤掉，以保证分析的正常进行。");
$section->submenu("测序数据质量分析",-icon => "./src/icon/2.1guolv.png");

# 2.1 过滤前后碱基质量图
my @fig_tree = ("./1.qc/${sample}.raw.png", "./1.qc/${sample}.cln.png");
my @labels = qw(过滤前 过滤后);
Image($section, \@fig_tree, \@labels, "碱基质量图");

# 2.2 数据统计
$section->submenu("测序数据统计",-icon => "./src/icon/2.2stat.png");
if (-s  "${dir}/1.qc/${sample}.nanostat" ){
    $section->tsv2html(
        -file   => "${dir}/1.qc/${sample}.nanostat",
        -top => 100,
        -name   => "数据统计结果",
        -header => 1);
}else{
    $section->desc("$blank6 无数据统计结果");}


# 3 拼接序列
my $section;
$section = $report->section(id => "consensus");
$section->menu("拼接序列",-icon => "./src/icon/3.ass.png");

# 3.1 拼接基因组
$section->submenu("组装序列",-icon => "./src/icon/3.1dna.png");
FileLink($section, "$blank6  拼接序列请点击  <link=./4.cons/consensus.fa>");

# 3.2 参考基因组覆盖度和深度统计
$section->submenu("参考基因组覆盖度和深度统计",-icon => "./src/icon/3.2ass_qc.png");
$section->tsv2html(
    -file   => "${dir}/2.align/coverage.tsv",
    -top    => 100,
    -name   => "拼接覆盖度深度信息",
    -header => 1);

# 3.3 参考基因组覆盖图
$section->submenu("参考基因组覆盖图",-icon => "./src/icon/4.1length.png");
$section->img2html(
    -file => "./2.align/coverage.png",
    -name => "组装序列长度分布图");


# 4 变异统计
my $section;
$section = $report->section(id => "mutations");
$section->menu("变异统计表",-icon => "./src/icon/4.3shijunti.png");
$section->tsv2html(
    -file   => "${dir}/3.muts/muts.tsv",
    -top    => 100,
    -name   => "拼接覆盖度深度信息",
    -header => 1);



# 生成HTML报告
$report->write();