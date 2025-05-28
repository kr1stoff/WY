#!/usr/bin/env perl

use lib "/sdbb/share/lib/GDHR";
use strict;
use GDHR;
use Utils;
use File::Basename;

# 命令函参数
my $prog = basename($0);
my $usage = "Usage:\n\t$prog <workdir> <sample_name>\n";
die $usage if $#ARGV != 1;
my $workdir = $ARGV[0];
my $sample = $ARGV[1];

# 标题
my $report;
$report = GDHR->new(-outdir  => $workdir,
                    -pipe    => "新冠全基因组测序分析报告",
                    -nonlazy => 1);
# 项目概述
my $section;
$section = $report->section(id => "introduction");
$section->menu("项目概述", -icon => "source/src/icons/gaishu.png");
$section->desc("SARS-CoV-2是冠状病毒科的新发病原体，又称2019新型冠状病毒（2019-nCoV），具有高传染、高隐蔽性以及高变异性的特点。这种特点给毒株的分析以及溯源带来的极大的困难，通过普通的宏基因组测序难以得到准确的分型与溯源数据。\n本项目通过扩增子建库，拼接获得新冠病毒全长基因组，可进行病毒株系的分型。同时，进行基于全长基因组的多序列比对进化树与基于突变的 SNPs 进化树分析，得到溯源数据。");

# 信息分析流程
$section = $report->section(id => "analysis");
$section->menu("信息分析流程", -icon => "source/src/icons/liucheng.png");
$section->desc("生物信息分析流程图，如下图所示");
$section->img2html(
    -file => "source/src/sc2_analysis.png",
    -name => "生物信息分析流程图");

# 数据质控
my $section;
$section = $report->section(id => "dataQC");
$section->menu("数据质控", -icon => "source/src/icons/cexuzhikong.png");
FileLink($section, "全部原始数据质控结果请点击 <link=./1.qc>");
## 过滤前fastqc
$section->submenu("测序数据质量分析", -icon => "source/src/icons/guolv.png");
$section->desc("过滤前碱基质量值图，展示原始下机序列上每个位置碱基的平均质量值，结果如下图所示 ");
if(-e "$workdir/1.qc/$sample.2_before_per_base_quality.png"){ # 双端测序两个图
    my @name = [ "imgs2html2" ];
    $section->img2html2(
        -file1 => "1.qc/$sample.1_before_per_base_quality.png",
        -name1 => "双端测序序列1碱基质量图",
        -file2 => "1.qc/$sample.2_before_per_base_quality.png",
        -name2 => "双端测序序列2碱基质量图",
        -names => "img2html2");
}else{ # 单端测序一个图
    $section->img2html(
        -file => "1.qc/$sample.1_before_per_base_quality.png",
        -name => "单端测序序列碱基质量图");
};
## 过滤后fastqc
$section->desc("过滤后碱基质量值图，展示过滤后序列上每个位置碱基的平均质量值，结果如下图所示 ");
if(-e "$workdir/1.qc/$sample.2_after_per_base_quality.png"){ # 双端测序两个图
    my @name = [ "imgs2html2" ];
    $section->img2html2(
        -file1 => "1.qc/$sample.1_after_per_base_quality.png",
        -name1 => "双端测序序列1碱基质量图",
        -file2 => "1.qc/$sample.2_after_per_base_quality.png",
        -name2 => "双端测序序列2碱基质量图",
        -names => "img2html2");
}else{ # 单端测序一个图
    $section->img2html(
        -file => "1.qc/$sample.1_after_per_base_quality.png",
        -name => "单端测序序列碱基质量图");
};
## 数据过滤前后的对比结果
$section->submenu("数据统计", -icon => "source/src/icons/guolv.png");
$section->tsv2html( # 添加一个表格
    -file   => "$workdir/1.qc/$sample.basic.stat.txt",
    -top    => 100,
    -name   => "过滤指标统计",
    -header => 1);

# 结果分析
my $section;
$section = $report->section(id => "results");
$section->menu("结果分析", -icon => "source/src/icons/fenxi.png");
## 比对结果
$section->submenu("比对结果", -icon => "source/src/icons/bidui.png");
FileLink($section, "全部比对质控结果请点击 <link=./2.map>");
$section->desc("将过滤后 FASTQ 比对回参考基因组，对初次比对结果的 InDel 区域进行重比对矫正比对错误，得到最终比对结果。结果如下表所示：");
$section->tsv2html(
    -file   => "$workdir/2.map/$sample.bam_stats.txt",
    -top    => 100,
    -name   => "基因比对覆盖度与深度表",
    -header => 1);
$section->desc("结果说明<br>平均深度：基因组平均覆盖深度，覆盖度：基因组覆盖度，深度≥10x：覆盖深度达到10层的区域在基因组中的占比，均一性：基因组覆盖深度的一致性指数。");
$section->desc("全基因组覆盖图如下图所示 ");
my @name = [ "imgs2html2" ];

$section->img2html(
-file => "2.map/$sample.genome_coverage_depth.png",
-name => "全基因组覆盖度图");

$section->img2html(
-file => "2.map/$sample.genome_coverage_depth_ylim1000.png",
-name => "全基因组覆盖度图（纵坐标限制最大值为1000");

$section->desc("结果说明: 上图为全基因组覆盖度图，下图为纵坐标限制最大值为1000的全基因组覆盖图。横坐标为参考基因组位置区间，纵坐标为该区间的平均测序覆盖深度，红色虚线为均一性阈值线。");
## 新冠病毒分型结果
$section->submenu("新冠病毒分型结果", -icon => "source/src/icons/xinguan.png");
$section->desc("具体谱系分型结果如下表所示:");
$section->tsv2html(
    -file   => "$workdir/source/$sample.lineage.tsv",
    -top    => 100,
    -name   => "新冠谱系分型结果表",
    -header => 1);
$section->desc("附：本次分析的所有样本，新冠病毒分型结果汇总，如下表所示：");
$section->tsv2html(
    -file   => "$workdir/source/lineage_report_trans.xls",
    -top    => 100,
    -name   => "本批次分析新冠谱系分型结果表",
    -header => 1);
FileLink($section, "新冠谱系分型结果查看请点击 <link=./source/lineage_report_trans.xlsx>");
## 样本 SNP/InDel 变异结果
$section->submenu("样本 SNP/InDel 变异结果", -icon => "source/src/icons/bianyi.png");
$section->desc("本分析流程的变异检测，是基于贝叶斯单倍型的遗传多态性算法。结果如下表所示：");
$section->tsv2html(
    -file   => "$workdir/3.variant/$sample.trans.tsv",
    -top    => 100,
    -name   => "SNP/InDel 结果表",
    -header => 1);
FileLink($section, "变异文件查看请点击 <link=./3.variant/$sample.trans.xlsx>");
# ## 进化树分析
$section->submenu("进化树分析");
$section->desc("系统进化树展示具有共同祖先的各物种间进化关系的树，是一种亲缘分支分类方法。在树中，每个节点代表其各分支的最近共同祖先，节点间层级关系代表物种间的亲缘关系。");
### 单样本多序列比对进化树
$section->ssubmenu("单样本多序列比对进化树");
$section->desc("以新冠病毒作参考基因组，进行比对，获得各样本比对信息，使用比对信息构建一致性序列，达到完整基因组拼接的效果。\n单个样本与 WHO 的 VOC/VOI 的代表性序列进行多序列比对，使用最大似然法构建进化树。如下图所示：");
if (-e "$workdir/5.phylogenetic/rectangular.png"){
    $section->img2html(
    -file => "5.phylogenetic/rectangular.png",
    -name => "进化树图");
}else{
    $section->desc("由于样本间SNP差异较小，本次分析无基于SNPs方法构建进化树结果。");
};
### SNP进化树
$section->ssubmenu("基于SNPs方法构建进化树");
$section->desc("比对参考基因组获得各样本变异信息，多样本间共有 SNP 生成特有标记性序列，使用最大似然法构建进化树。");
$section->desc("使用 SNP-based 方法构建进化树相较于“基于多序列比对方法构建进化树”优势是在原始数据序列完整度不足时也可实现进化树的构建。");
if (-e "$workdir/source/SNPTree/rectangular.png"){
    $section->img2html(
    -file => "source/SNPTree/rectangular.png",
    -name => "进化树图");
}else{
    $section->desc("由于样本间SNP差异较小，本次分析无基于SNPs方法构建进化树结果。");
};
### MSA进化树
$section->ssubmenu("基于多序列比对方法构建进化树");
$section->desc("比对参考基因组获得各样本变异信息，使用比对信息构建一致性序列，达到完整基因组拼接的效果。多个样本的完整基因组进行多序列比对，最终使用最大似然法构建进化树（ML Phylogenetic tree）。 使用基于多序列比对方法构建进化树方法相较于“基于SNPs方法构建进化树”有更高的准确性。");
if (-e "$workdir/source/MSATree/rectangular.png"){
    $section->img2html(
    -file => "source/MSATree/rectangular.png",
    -name => "进化树图");
}else{
    $section->desc("由于样本间差异较小，本次分析无基于多序列比对方法构建进化结果。");
};
## 新冠全球库
if (-e "$workdir/6.global_database"){ # 可能没有选
    $section->submenu("新冠全球库分析", -icon => "source/src/icons/quanqiu.png");
    $section->desc("使用拼接好的完整基因组与本地新冠全球库进行比对，选取相似度最高的部分序列，进行溯源进化树分析。目前系统本地中收录病毒毒株数超五百万条，并持续更新。");
    if (-e "$workdir/6.global_database/rectangular.png"){
        $section->img2html(
        -file => "6.global_database/rectangular.png",
        -name => "进化树图");
    }else{
        $section->desc("由于样本间差异较小，本次分析无基于多序列比对方法构建进化结果。");
    }
};

$section->menu("检测局限性", -icon => "source/src/icons/gaishu.png");
$section->desc("1.检测结果仅供科研使用，如有其它疑义请在七个工作日内与我们联系。");
$section->desc("2.样本采集，保存，运输等环节的不当操作可能影响本产品检测性能。");
$section->desc("3.检测方将严格依法保护患者隐私与检测结果。<br><br>");
$section->desc("<font color='blue'>我们深深明了，每一份珍贵样本的背后，都有一位亟待救治的病患</font><br>");

$section->desc("<br>总 部： 中国·广州 黄埔区华南新材料创新园G10栋<br>公司电话：400-809-9973 | 网址：visionmedicals.com | 邮箱：IDseq\@visionmedicals.com<br>");

# 生成HTML报告
$report->write();
