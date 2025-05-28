#!/usr/bin/env perl
my $app;
my $GDHRPATH;
BEGIN {
    use FindBin qw($Bin $RealBin);
    use YAML qw(LoadFile);
    my $f_software = "$RealBin/../config/software.yml";
    $app = LoadFile($f_software);
    $GDHRPATH = $app->{"GDHR"};
}

use strict;
use Getopt::Long;
use Cwd qw(realpath);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib $GDHRPATH;
use GDHR;
use Utils;

#### param
my %OPTS = (
    nonlazy => 1,
    company => "微远基因",
    outdir  => "./"
);
GetOptions(\%OPTS,
    'outdir:s',
    'nonlazy:s',
    'company:s',
    'help'
);
&usage if ($OPTS{help});

## Some global variable
$OPTS{outdir} = realpath($OPTS{outdir});
our $dir = ".";
our $fulldir = $OPTS{outdir};
my @segments = qw/PB2 PB1 PA HA NP NA MP NS/;

#### Main
my $report = GDHR->new(
    -outdir  => $OPTS{outdir},
    -company => "微远基因",
    -pipe    => "流感扩增子测序分析报告",
    nonlazy  => $OPTS{nonlazy}
);

&qc($report);
&assemble($report);
&coverage($report);
our $f_assem_stat = "$fulldir/02.Assemble/assemble.stat.tsv";
our $is_ha = &is_ha($f_assem_stat);
our $is_na = &is_na($f_assem_stat);
if ($is_ha + $is_na > 0) {
    &blast($report);
}
&phylotree($report);
if ($is_ha + $is_na > 0) {
    &snp($report);
}
$report->write();
`cp $Bin/../etc/image/* $fulldir/src/image`;
$report->pack(-format => "zip");

#### Some Function
sub qc() {
    my $report = shift;

    my $dir0 = $dir;
    my $fulldir0 = $fulldir;
    my $dir = "$dir0/01.QC";
    my $fulldir = "$fulldir0/01.QC";

    my $section = $report->section(id => "qc");
    $section->menu("数据质控", -icon => "./src/image/qc.png");
    my $desc_qc = <<CONT;
测序数据下机后，我们会对其进行质控，筛选出高质量的 reads 进行后续分析，质控后数据碱基质量分布图如下：
CONT
    $section->desc($desc_qc);
    if (-e "$fulldir/qc_2.png") {
        $section->img2html2(
            -file1 => "$dir/qc_1.png",
            -name1 => "Read1",
            -file2 => "$dir/qc_2.png",
            -name2 => "Read2");
    }
    else {
        $section->img2html(
            -file => "$dir/qc.png",
            -name => "After QC");
    }
    $section->desc("过滤前后数据统计如下");
    $section->tsv2html(
        -file   => "$fulldir/qc.stat.tsv",
        -name   => "过滤前后数据统计",
        -header => 1);
}

sub assemble() {
    my $report = shift;

    my $dir0 = $dir;
    my $fulldir0 = $fulldir;
    my $dir = "$dir0/02.Assemble";
    my $fulldir = "$fulldir0/02.Assemble";

    my $section = $report->section(id => "assemble");
    $section->menu("序列组装", -icon => "./src/image/assemble.png");
    $section->desc("流感片段长度统计");
    $section->tsv2html(
        -file   => "$fulldir/assemble.stat.tsv",
        -name   => "流感片段长度统计",
        -header => 1);
    my @sequence = map {"$dir/${_}.fasta"} @segments;
    my @desc_sequence = map {"${_}"} @segments;
    $section->desc(" 序列文件下载");
    $section->files2list(
        -files     => \@sequence,
        -desc      => \@desc_sequence,
        -short_dir => \@desc_sequence);
}

sub coverage() {
    my $report = shift;

    my $dir0 = $dir;
    my $fulldir0 = $fulldir;
    my $dir = "$dir0/07.Coverage";
    my $fulldir = "$fulldir0/07.Coverage";

    my $section = $report->section(id => "coverage");
    $section->menu("覆盖度统计", -icon => "./src/image/coverage.png");

    $section->tsv2html(
        -file   => "$fulldir/coverage.tsv",
        -name   => "覆盖度统计",
        -header => 1);

    my @imgs1 = map {"$dir/${_}.genome_coverage_depth.png"} @segments;
    my @imgs2 = map {"$dir/${_}.genome_coverage_depth_ylim1000.png"} @segments;
    my @names1 = map {"${_}覆盖度图"} @segments;
    my @names2 = map {"${_}覆盖度图(纵坐标限制最大值为1000)"} @segments;
    $section->imgs2html2(
        -files1 => \@imgs1,
        -files2 => \@imgs2,
        -desc1  => \@names1,
        -desc2  => \@names2,
        -names  => \@segments);
}

sub blast() {
    my $report = shift;

    my $dir0 = $dir;
    my $fulldir0 = $fulldir;
    my $dir = "$dir0/03.Blast";
    my $fulldir = "$fulldir0/03.Blast";

    my $section = $report->section(id => "blast");
    $section->menu("序列注释", -icon => "./src/image/blast.png");
    if (&is_ha($f_assem_stat)) {
        $section->tsv2html(
            -file   => "$fulldir/HA.annot.xls",
            -name   => "HA基因序列注释结果",
            -header => 1);
    }
    if (&is_na($f_assem_stat)) {
        $section->tsv2html(
            -file   => "$fulldir/NA.annot.xls",
            -name   => "NA基因序列注释结果",
            -header => 1);
    }
}

sub phylotree() {
    my $report = shift;

    my $dir0 = $dir;
    my $fulldir0 = $fulldir;
    my $dir = "$dir0/04.PhyloTree";
    my $fulldir = "$fulldir0/04.PhyloTree";

    my $section = $report->section(id => "phylotree");
    $section->menu("进化树分析", -icon => "./src/image/phylotree.png");
    $section->desc("系统进化树展示具有共同祖先的各物种间进化关系的树，是一种亲缘分支分类方法。");
    my @imgs = map {"$dir/${_}.png"} @segments;
    $section->imgs2html(
        -files => \@imgs,
        -names => \@segments,
        -name  => "流感片段进化树分析 ");
}

sub snp() {
    my $report = shift;

    my $dir0 = $dir;
    my $fulldir0 = $fulldir;
    my $dir = "$dir0/06.SNP";
    my $fulldir = "$fulldir0/06.SNP";

    my $section = $report->section(id => "snp");
    $section->menu("突变分析", -icon => "./src/image/snp.png");

    if (&is_ha($f_assem_stat)) {
        $section->tsv2html(
            -file   => "$fulldir/HA.variants.txt",
            -name   => "HA基因突变信息",
            -header => 1);

    }
    if (&is_na($f_assem_stat)) {
        $section->tsv2html(
            -file   => "$fulldir/NA.variants.txt",
            -name   => "NA基因突变信息",
            -header => 1);
    }

}

sub is_ha() {
    my $f_assem_stat = shift;

    open(INPUT, "<$f_assem_stat") or die "$f_assem_stat 文件无法打开, $!";
    my $flag = 0;
    while (<INPUT>) {
        if ($flag == 1) {
            my @arr = split("\t", $_);
            if ($arr[3] == 0) {
                return (0);
            }
            else {
                return (1);
            }
        }
        else {
            $flag += 1;
        }
    }
    return (1)
}

sub is_na() {
    my $f_assem_stat = shift;

    open(INPUT, "<", $f_assem_stat) or die "$f_assem_stat 文件无法打开, $!";
    my $flag = 0;
    while (<INPUT>) {
        if ($flag == 1) {
            my @arr = split("\t", $_);
            if ($arr[5] == 0) {
                return (0);
            }
            else {
                return (1);
            }
        }
        else {
            $flag += 1;
        }
    }
    return (1)
}

sub usage {
    my $help = <<HELP;

Useage:
         perl $0 -outdir Pipe

Options: -outdir  STR     the out put dir [required]
         -nonlazy INT     lazy or non-lazy [default 1]
         -company STR     The company name [default 微远基因]
         -help|h          output this information

HELP
    print $help;
    exit 0;
}



