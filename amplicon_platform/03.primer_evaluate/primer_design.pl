#!/usr/bin/perl -w 
#########################################################
# Author: stone
# Created Time : 2023.05.25 10:00:00 AM CST
# File Name: primer_design.pl
# Version: v1.0
#########################################################

use strict;
use warnings;
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use FindBin qw($Bin);
# use lib "$Bin/../lib";
use lib "$Bin/lib";
use GACP qw(parse_config);


# 帮助文档
=head1 Description

=encoding utf-8

    此软件用于多重PCR引物设计。

=head1 Usage

    perl $0 -regions <bed file> -type <qPCR|tNGS|Both>

=head1 Parameters

    optional arguments:
    -regions                   bed文件（染色体\t起始\t终止\t物种拉丁名-xxx\t分管信息\tB\t整体设计[W]或非整体设计[NotW]）
    -draft                     草稿引物
    -ref                       参考基因组
    -minampllen                扩增子最小长度，默认175
    -maxampllen                扩增子最大长度，默认250
    -optampllen                扩增子最优长度，默认200
    -minprimerlen              引物最小长度，默认18
    -maxprimerlen              引物最大长度，默认28
    -optprimerlen              引物最优长度，默认23
    -minprimermelt             引物最小退火温度，默认63
    -maxprimermelt             引物最大退火温度，默认67
    -optprimermelt             引物最优退火温度，默认64
    -minprimergc               引物最小GC含量，默认30
    -maxprimergc               引物最大GC含量，默认70
    -optprimergc               引物最优GC含量，默认50
    -minprimerendgc            引物3'端末尾5个碱基含G或C的最小数量，默认0
    -maxprimerendgc            引物3'端末尾5个碱基含G或C的最大数量，默认5
    -optprimerendgc            引物3'端末尾5个碱基含G或C的最优数量，默认2
    -maxprimerpolyn            引物中最大可接受的相同连续碱基的数量，默认5
    -maxoverlap                扩增子之间的overlap碱基最大数量，默认300
    -primernum1                扩增子的引物数量，默认20
    -autoadjust                引物设计不成功，自动调整参数（T调整，F不调整,默认F）
    -returnvariantsnum         输出引物组合的数目，默认5
    -blast                     引物比对，默认TRUE
    -snps                      引物设计排除含有SNP的点（可选）
    -dbsnp                     突变信息VCF文件
    -nucs                      引物设计排除3'端多少数量的碱基含有SNP的点，默认引物全部碱基
    -th                        线程数，默认100
    -run                       输出文件的前缀名
    -skip                      跳过不能用默认参数设计引物的区域（可选）
    -type                      设计引物的类型(必填之一：tNGS、qPCR、Both), 默认tNGS
    -outdir                    输出结果的目录
    -genome_dir                输入数据的目录
    -help                      显示帮助信息

    #引物长度：最优：18-28bp，次优：17-30bp；最优长度：23bp；GC含量范围：最优30-70%，次优20-80%，最优GC：50%，次优GC：45%；退火温度：63-67℃，次优退火温度：62-67℃，最优：64℃；次优：65℃#

=cut

my ($type,$outdir,$genome_dir);
my ($regions,$draft,$ref,$minampllen,$maxampllen,$optampllen,$minprimerlen,$maxprimerlen,$optprimerlen,$minprimermelt,$maxprimermelt,$optprimermelt);
my ($minprimergc,$maxprimergc,$optprimergc,$minprimerendgc,$maxprimerendgc,$optprimerendgc,$maxprimerpolyn,$maxoverlap,$primernum1,$autoadjust);
my ($returnvariantsnum,$blast,$snps,$dbsnp,$nucs,$th,$run,$skip,$help);

GetOptions(
    "regions=s"                => \$regions,
    "draft:s"                  => \$draft,
    "ref:s"                    => \$ref,
    "minampllen:s"             => \$minampllen,
    "maxampllen:s"             => \$maxampllen,
    "optampllen:s"             => \$optampllen,
    "minprimerlen:s"           => \$minprimerlen,
    "maxprimerlen:s"           => \$maxprimerlen,
    "optprimerlen:s"           => \$optprimerlen,
    "minprimermelt:s"          => \$minprimermelt,
    "maxprimermelt:s"          => \$maxprimermelt,
    "optprimermelt:s"          => \$optprimermelt,
    "minprimergc:s"            => \$minprimergc,
    "maxprimergc:s"            => \$maxprimergc,
    "optprimergc:s"            => \$optprimergc,
    "minprimerendgc:s"         => \$minprimerendgc,
    "maxprimerendgc:s"         => \$maxprimerendgc,
    "optprimerendgc:s"         => \$optprimerendgc,
    "maxprimerpolyn:s"         => \$maxprimerpolyn,
    "maxoverlap:s"             => \$maxoverlap,
    "primernum1:s"             => \$primernum1,
    "autoadjust:s"             => \$autoadjust,
    "returnvariantsnum:s"      => \$returnvariantsnum,
    "blast:s"                  => \$blast,
    "snps:s"                   => \$snps,
    "dbsnp:s"                  => \$dbsnp,
    "nucs:s"                   => \$nucs,
    "th:s"                     => \$th,
    "run:s"                    => \$run,
    "skip:s"                   => \$skip,
    "type=s"                   => \$type,
    "outdir:s"                 => \$outdir,
    "genome_dir=s"             => \$genome_dir,
    "h|help:s"                 => \$help
    );

die `pod2text $0` if ((!$regions)  or ($help));

#默认参数
# my $config_file = "$Bin/../config.txt";
my $config_file = "$Bin/config.txt";
my $name2taxid_info = parse_config($config_file,"name2taxid_info");
my $all_genome_dir = parse_config($config_file,"all_genome_dir");
my $csvtk = parse_config($config_file,"csvtk");


$minampllen ||= 175;
$maxampllen ||= 250;
$optampllen ||= 200;
$minprimerlen ||= 18;
$maxprimerlen ||= 28;
$optprimerlen ||= 23;
$minprimermelt ||= 63;
$maxprimermelt ||= 67;
$optprimermelt ||= 64;
$minprimergc ||= 30;
$maxprimergc ||= 70;
$optprimergc ||= 50;
$minprimerendgc ||= 0;
$maxprimerendgc ||= 5;
$optprimerendgc ||= 2;
$maxprimerpolyn ||= 5;
$maxoverlap ||= 300;
$primernum1 ||= 10;
$returnvariantsnum ||= 5;
$th ||= 100;
$autoadjust ||= "F";
$type ||= "tNGS";

$outdir ||= "./primergo";
print "mdkir -p $outdir\n";
`mkdir -p $outdir`;
$outdir = abs_path($outdir);

$regions = abs_path($regions);


my %species2taxid;
open IN,$name2taxid_info or die $!;  # name to taxid info file
while(<IN>){
    chomp;
    next if /#/;
    my @line = split /\t/,$_;
    $species2taxid{$line[0]} = $line[1];
}
close IN;

#拷贝regions到结果目录
my %species;
my ($taxid,$spe_name);
open IN,"<$regions" or die $!;
while(<IN>){
    chomp;
    my @line = split/\t/,$_;
    my @name_tmps = split/-/,$line[3];
    $spe_name = $name_tmps[0];
    if ($species2taxid{$spe_name}){
        $taxid = $species2taxid{$spe_name};
    }
    else{
        die "bed文件中,物种名不对,或者物种名没有对应的taxid\n";
    }
}
close IN;

print "Species: $spe_name\n";
print "Taxid: $taxid\n";

$genome_dir ||= "$all_genome_dir/$taxid";
$ref ||= "$genome_dir/ref.fna";

`cp $regions $outdir/2tNGS.bed`;

if($type eq "Both"){
    my $s = `date "+%Y-%m-%d %H:%M:%S"`;
    print "\n开始引物设计:\n$s\n";
    if($autoadjust eq "T"){
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py -ref $ref -region $outdir/2tNGS.bed -th 100 -primernum1 $primernum1 -blast --embedded-amplification -minampllen 70 -maxampllen 110 -optampllen 100 -optextampllen 200 -maxextampllen 250 -returnvariantsnum 5 -skip -maxoverlap 300 -maxprimerpolyn $maxprimerpolyn -autoadjust`;
    }else{
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py -ref $ref -region $outdir/2tNGS.bed -th $th -primernum1 $primernum1 -blast --embedded-amplification -minampllen 70 -maxampllen 110 -optampllen 100 -optextampllen 200 -maxextampllen 250 -minprimerlen $minprimerlen -maxprimerlen $maxprimerlen -optprimerlen $optprimerlen -minprimermelt $minprimermelt -maxprimermelt $maxprimermelt -optprimermelt $optprimermelt -minprimergc $minprimergc -maxprimergc $maxprimergc -optprimergc $optprimergc -minprimerendgc $minprimerendgc -maxprimerendgc $maxprimerendgc -maxprimerpolyn $maxprimerpolyn  -maxoverlap $maxoverlap -returnvariantsnum $returnvariantsnum -skip -autoadjust`;
    }
    print "将excel转成txt\n";
    `$csvtk xlsx2csv -t $outdir/2tNGS_primers_combination_1_info.xls > $outdir/2tNGS_primers_combination_1_info.txt`;
    my $e = `date "+%Y-%m-%d %H:%M:%S"`;
    print "引物设计结束:\n$e\n";
}elsif($type eq "qPCR"){
    my $s = `date "+%Y-%m-%d %H:%M:%S"`;
    print "\n开始引物设计:\n$s\n";
    if($autoadjust eq "T"){
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py -ref $ref -region $outdir/2tNGS.bed -th $th -primernum1 $primernum1 -blast -minampllen 70 -maxampllen 110 -optampllen 100 -minprimerlen $minprimerlen -maxprimerlen $maxprimerlen -optprimerlen $optprimerlen -minprimermelt $minprimermelt -maxprimermelt $maxprimermelt -optprimermelt $optprimermelt -minprimergc $minprimergc -maxprimergc $maxprimergc -optprimergc $optprimergc -minprimerendgc $minprimerendgc -maxprimerendgc $maxprimerendgc -maxprimerpolyn $maxprimerpolyn -maxoverlap $maxoverlap -returnvariantsnum $returnvariantsnum -skip -autoadjust`;
    }else{
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py -ref $ref -region $outdir/2tNGS.bed -th $th -primernum1 $primernum1 -blast -minampllen 70 -maxampllen 110 -optampllen 100 -minprimerlen 15 -maxprimerlen 25 -optprimerlen 20 -minprimermelt 58 -maxprimermelt 62 -optprimermelt 60 -minprimergc $minprimergc -maxprimergc $maxprimergc -optprimergc $optprimergc -minprimerendgc $minprimerendgc -maxprimerendgc $maxprimerendgc -maxprimerpolyn $maxprimerpolyn  -maxoverlap $maxoverlap -returnvariantsnum $returnvariantsnum -skip`;
    }   
    print "将excel转成txt\n";
    `$csvtk xlsx2csv -t $outdir/2tNGS_primers_combination_1_info.xls > $outdir/2tNGS_primers_combination_1_info.txt`;
    my $e = `date "+%Y-%m-%d %H:%M:%S"`;
    print "引物设计结束:\n$e\n";
}elsif($type eq "tNGS"){
    my $s = `date "+%Y-%m-%d %H:%M:%S"`;
    print "\n开始引物设计:\n$s\n";
    if($autoadjust eq "T"){
        print "/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py  -ref $ref -region $outdir/2tNGS.bed -th 100 -primernum1 $primernum1 -blast -minampllen 175 -maxampllen 250 -optampllen 200 -returnvariantsnum 5 -skip -maxoverlap 300 -autoadjust -maxprimerpolyn $maxprimerpolyn\n";
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py  -ref $ref -region $outdir/2tNGS.bed -th 100 -primernum1 $primernum1 -blast -minampllen 175 -maxampllen 250 -optampllen 200 -returnvariantsnum 5 -skip -maxoverlap 300 -autoadjust -maxprimerpolyn $maxprimerpolyn`;
    }else{
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py -ref $ref -region $outdir/2tNGS.bed -th $th -primernum1 $primernum1 -blast -minampllen $minampllen -maxampllen $maxampllen -optampllen $optampllen -minprimerlen $minprimerlen -maxprimerlen $maxprimerlen -optprimerlen $optprimerlen -minprimermelt $minprimermelt -maxprimermelt $maxprimermelt -optprimermelt $optprimermelt -minprimergc $minprimergc -maxprimergc $maxprimergc -optprimergc $optprimergc -minprimerendgc $minprimerendgc -maxprimerendgc $maxprimerendgc -maxprimerpolyn $maxprimerpolyn  -maxoverlap $maxoverlap -returnvariantsnum $returnvariantsnum -skip`;
    }
    print "将excel转成txt\n";
    `$csvtk xlsx2csv -t $outdir/2tNGS_primers_combination_1_info.xls > $outdir/2tNGS_primers_combination_1_info.txt`;
    my $e = `date "+%Y-%m-%d %H:%M:%S"`;
    print "引物设计结束:\n$e\n";
}
elsif($type eq "seq"){
    my $s = `date "+%Y-%m-%d %H:%M:%S"`;
    print "\n开始引物设计:\n$s\n";
    if($autoadjust eq "T"){
        print "/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py  -ref $ref -region $outdir/2tNGS.bed -th 100 -primernum1 $primernum1 -blast -minampllen 280 -maxampllen 500 -optampllen 300 -returnvariantsnum 5 -skip -maxoverlap 300 --min-primer-length 17 --max-primer-length 28 --optimal-primer-length 23 --min-primer-gc 40 --max-primer-gc 60 --optimal-primer-gc 50 --min-primer-melting-temp 63 --max-primer-melting-temp 67 --optimal-primer-melting-temp 64 -autoadjust -maxprimerpolyn $maxprimerpolyn\n";
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py  -ref $ref -region $outdir/2tNGS.bed -th 100 -primernum1 $primernum1 -blast -minampllen 280 -maxampllen 500 -optampllen 300 -returnvariantsnum 5 -skip -maxoverlap 300 --min-primer-length 17 --max-primer-length 28 --optimal-primer-length 23 --min-primer-gc 40 --max-primer-gc 60 --optimal-primer-gc 50 --min-primer-melting-temp 63 --max-primer-melting-temp 67 --optimal-primer-melting-temp 64 -autoadjust -maxprimerpolyn $maxprimerpolyn`;
    }else{
        `/home/yehui/software/python3/Python-3.8.15/bin/python3.8 /sdbb/bioinfor/mengxf/Software/github/NGS-PrimerPlex-v1.3.4/NGS_primerplex.py -ref $ref -region $outdir/2tNGS.bed -th $th -primernum1 $primernum1 -blast -minampllen $minampllen -maxampllen $maxampllen -optampllen $optampllen -minprimerlen $minprimerlen -maxprimerlen $maxprimerlen -optprimerlen $optprimerlen -minprimermelt $minprimermelt -maxprimermelt $maxprimermelt -optprimermelt $optprimermelt -minprimergc $minprimergc -maxprimergc $maxprimergc -optprimergc $optprimergc -minprimerendgc $minprimerendgc -maxprimerendgc $maxprimerendgc -maxprimerpolyn $maxprimerpolyn  -maxoverlap $maxoverlap -returnvariantsnum $returnvariantsnum -skip`;
    }
    print "将excel转成txt\n";
    `$csvtk xlsx2csv -t $outdir/2tNGS_primers_combination_1_info.xls > $outdir/2tNGS_primers_combination_1_info.txt`;
    my $e = `date "+%Y-%m-%d %H:%M:%S"`;
    print "引物设计结束:\n$e\n";
}else{
    print "请检查设计引物的类型！！！\n";
}

