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
$section->desc("$blank6 通过分离培养获得病毒株，进行建库测序，获得较为纯净的毒株序列。本项目对毒株（样本编号为$sample）的Nanopore测序数据进行分析，通过从头组装最终得到高质量的基因组序列，并对该序列进行基因预测及基因注释。");


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
$section->submenu("数据统计",-icon => "./src/icon/2.2stat.png");
if (-s  "${dir}/1.qc/${sample}.nanostat" ){
    $section->tsv2html(
        -file   => "${dir}/1.qc/${sample}.nanostat",
        -top => 100,
        -name   => "数据统计结果",
        -header => 1);
}else{
    $section->desc("$blank6 无数据统计结果");}


# 3 组装基因组
$section = $report->section(id => "assemble");
$section->menu("组装基因组",-icon => "./src/icon/3.ass.png");
if (-e "${dir}/${sample}.fa"){

    # 3.1 组装序列
    FileLink($section, "$blank6  全部组装质控结果请点击  <link=./2.asm>");
    $section->submenu("组装序列",-icon => "./src/icon/3.1dna.png");
    FileLink($section, "$blank6  组装序列请点击  <link=./${sample}.fa>");
    
    # 3.2.1 基本组装指标
    $section->submenu("组装质控",-icon => "./src/icon/3.2ass_qc.png");
    $section->ssubmenu("基本组装指标");
    $section->desc("$blank6 组装完成后，可得到一系列用于表征组装质量的参数，包括scaffolds数目，scaffolds总长度、N50等。一般来说，scaffolds数目越少，N50、N90值越高，表示组装效果越好。");
    if (-e "${dir}/2.asm/${sample}_fa.stat.txt"){
    $section->tsv2html(
        -file   => "${dir}/2.asm/${sample}_fa.stat.txt",
        -top    => 100,
        -name   => "${sample}组装质控信息",
        -header => 1);
    }else{
        $section->desc("$blank6 无组装质控信息");}
    $section->desc("$blank6 注：N50：将组装序列按照长度，从长到短进行排序，将序列长度依次相加，当累计长度达到组装序列总长度50%时，当前累加的序列长度即为 N50，累加的contig的个数即为L50，同理可得N75与N90，L75与L90。");

    $section->submenu("组装质控",-icon => "./src/icon/3.2ass_qc.png");

    # 3.2.2 组装序列长度分布图
    $section->ssubmenu("组装基因组长度分布图");
    if (-e "${dir}/2.asm/${sample}.length.png"){
    $section->img2html(
        -file => "./2.asm/${sample}.length.png",
        -name => "组装序列长度分布图");
    }else{
        $section->desc("$blank6 无组装序列长度分布图");}

    # 3.2.3 一致性
    $section->ssubmenu("组装基因组一致性");
    $section->desc("$blank6 将抽取reads比对回组装序列上，统计抽取reads对序列的覆盖情况，从而评估组装结果与测序数据的一致性。较高的覆盖率（95%以上）认为组装结果和reads有比较好的一致性。");
    if (-s  "${dir}/2.asm/consistency.txt"){
        $section->tsv2html(
            -file   => "${dir}/2.asm/consistency.txt",
            -top    => 100,
            -name   => "组装一致性统计结果",
            -header => 1);
        FileLink($section, "$blank6  每条scaffolds序列的组装一致性统计结果请点击  <link=./2.asm/consistency.all.txt>");
    }else{
        $section->desc("$blank6 无组装一致性统计结果");}

    # 3.2.4 GC分布
    $section->ssubmenu("组装基因深度分布");
    $section->desc("$blank6 通过用滑窗法切割组装基因组，来评估测序是否有明显的GC偏向，也可判断样品是否存在污染等情况。如果多数的点集中分布在一个比较窄的范围内，则表明样本不存在污染；如果分布在多个区域，则表明样本中可能存在其它物种的污染。当存在污染时，可根据点分布的集中区判断污染程度，例如根据GC分布大致判断污染物种的数量，或者根据污染序列的测序reads覆盖深度来大致推测污染的比例等");
    $section->img2html(
        -file => "./2.asm/depth_base.stat.depth_GC.png",
        -name => "组装基因深度分布");
    $section->desc("注： 散点图，横坐标为GC含量，纵坐标为reads覆盖深度；两侧散点图，分别为GC含量、测序深度的滑窗频数分布。通过GC_depth分布图可以看出测序是否有明显的GC偏向。");

    # 3.2.5 checkv
    $section->ssubmenu("组装基因组完整度/污染度评估");
    $section->desc("$blank6 一个谱系(lineage)的基因组会有一些共享的基因，称之为lineage-specific marker set。");
    $section->desc("$blank6 评估一个基因组的完整度，首先判断这个基因组的lineage，通过检索该lineage的marker set，如果marker set的基因全部找到了，则认为这个基因组完整度100%，若检测到多个marker，则有可能是杂合或者污染。");
    $section->desc("$blank6 本次分析对病毒基因组组装的完整度、杂合度质量评估");

    if (-s  "${dir}/2.asm/quality_summary.tsv"){
        $section->tsv2html(
            -file   => "${dir}/2.asm/quality_summary.tsv",
            -top    => 100,
            -name   => "组装完整度和准确度统计",
            -header => 1);
    }else{
        $section->desc("$blank6 无组装完整度和准确度统计结果");}
}


# 4 基因预测
$section = $report->section(id => "predict");
$section->menu("基因预测",-icon => "./src/icon/4.pre.png");
FileLink($section, "$blank6 全部基因预测结果请点击  <link=./3.pred>");

if (-e "${dir}/3.pred/predict.length.png"){
    # 4.1 基因长度分布图
    $section->submenu("预测基因长度分布",-icon => "./src/icon/4.1length.png");
    $section->img2html(
        -file => "./3.pred/predict.length.png",
        -name => "预测基因的长度分布图");
    $section->desc("$blank6 横坐标为基因长度范围（bp），纵坐标为对应长度范围的预测的基因数量");

    # 4.2 预测基因类型统计
    $section->submenu("预测基因分类统计",-icon => "./src/icon/4.2class.png");
    if (-s "${dir}/3.pred/predict_kind.txt"){
        $section->tsv2html(
            -file   => "${dir}/3.pred/predict_kind.txt",
            -top    => 100,
            -name   => "预测基因分类统计表",
            -header => 1);
    }else{
        $section->desc("$blank6 无预测基因分类统计表结果");}
}


# 5 基因注释
$section = $report->section(id => "anno");
$section->menu("基因注释",-icon => "./src/icon/5.anno.png");
FileLink($section, "$blank6  全部基因注释结果请点击  <link=./4.anno>");
# KEGG注释
$section->submenu("KEGG注释",-icon => "./src/icon/5.6kegg.png");
$section->desc("$blank6 基因会参与人体的各个通路。KEGG就是基于人体通路而形成的其中一个数据库，由日本京都大学生物信息学中心的Kanehisa实验室于1995年建立,数据库大致分为系统信息、基因组信息和化学信息三大类。");

# KEGG注释分类图
if (-s "${dir}/4.anno/KEGG_anno.stat.tsv"){
    $section->tsv2html(
        -file   => "${dir}/4.anno/KEGG_anno.stat.tsv",
        -top    => 100,
        -name   => "KEGG注释结果",
        -header => 1);
    $section->img2html(
        -file => "./4.anno/KEGG_anno_stats.png",
        -name => "KEGG注释分类图");
    FileLink($section, "$blank6  详细KEGG注释分类结果点击  <link=./4.anno/KEGG_anno.stat.tsv>");
}else{
    $section->desc("$blank6 无KEGG注释结果");}

# swiss-prot注释
$section->submenu("swiss-prot注释",-icon => "./src/icon/5.7swiss.png");
$section->desc("$blank6 SWISS-PROT是经过注释的蛋白质序列数据库，由欧洲生物信息学研究所(EBI)维护。数据库由蛋白质序列条目构成，每个条目包含蛋白质序列、引用文献信息、分类学信息、注释等，注释中包括蛋白质的功能、转录后修饰、特殊位点和区域、二级结构、四级结构、与其它序列的相似性、序列残缺与疾病的关系、序列变异体和冲突等信息。");
if (-s "${dir}/4.anno/swissprot_result.tsv"){
    $section->tsv2html(
        -file   => "${dir}/4.anno/swissprot_result.tsv",
        -top    => 100,
        -name   => "swiss-prot注释",
        -header => 1);
    FileLink($section, "$blank6 全部swiss-prot注释结果点击  <link=./4.anno/swissprot_result.tsv>");
}else{
    $section->desc("$blank6 无swiss-prot注释结果");}


# 生成HTML报告
$report->write();