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
    'species:s',
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
my $blank6 = blank(6);


#### Main
my $report = GDHR->new(
    -outdir  => $OPTS{outdir},
    -company => "微远基因",
    -pipe    => "$OPTS{species}分型",
    nonlazy  => $OPTS{nonlazy}
);

&mlst($report);
&cgmlst($report);
&wgmlst($report);
&cgsnp($report);
&wgsnp($report);
&typing($report);
&mlva($report);
&cansnp($report);
&spa($report);
`cp -rf $Bin/../etc/image/* $fulldir/src/image`;
$report->write();


#### Some Function
sub mlst {
    # MLST 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $file = "$fulldir/mlst.tsv";

    if (-e $file) {
        my $section = $report->section(id => "mlst");
        $section->menu("MLST 分型", -icon => "./src/image/mlst.png");
        my $desc_analysis = <<CONT;
$blank6
CONT
        $section->desc($desc_analysis);
        $section->tsv2html(
            -file   => $file,
            -name   => "MLST 分型结果",
            -header => 1);
    }
}

sub cgmlst {
    # cgMLST 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $f_info = "$fulldir/cgMLST.match.info.tsv";
    my $f_stat = "$fulldir/cgMLST.stat.tsv";

    if (-e $f_info) {
        my $section = $report->section(id => "cgmlst");
        $section->menu("cgMLST 分型", -icon => "./src/image/cgmlst.png");
        my $desc_analysis = <<CONT;
完全一致的基因如下表所示：
CONT
        $section->desc($desc_analysis);
        $section->tsv2html(
            -file   => $f_info,
            -name   => "完全匹配基因信息(Top10)",
            -header => 1,
            -top    => 10);

        if (-e $f_stat) {
            `sed -i '1i\\NA\\tNA' $f_stat`;
            $section->desc("cgMLST匹配信息如下：");
            $section->tsv2html(
                -file   => $f_stat,
                -name   => "cgMLST 结果统计",
                -header => 0);
            `sed -i '1d' $f_stat`;
        }

        FileLink($section, "canMLST详细结果 <link=./cgMLST.match.info.tsv>");

    }

}

sub wgmlst {
    # wgMLST 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $f_info = "$fulldir/wgMLST.match.info.tsv";
    my $f_stat = "$fulldir/wgMLST.stat.tsv";

    if (-e $f_info) {
        my $section = $report->section(id => "wgmlst");
        $section->menu("wgMLST 分型", -icon => "./src/image/wgmlst.png");
        my $desc_analysis = <<CONT;
完全一致的基因如下表所示：
CONT
        $section->desc($desc_analysis);
        $section->tsv2html(
            -file   => $f_info,
            -name   => "完全匹配基因信息(Top10)",
            -header => 1,
            -top    => 10);

        if (-e $f_stat) {
            `sed -i '1i\\NA\\tNA' $f_stat`;
            $section->desc("wgMLST匹配信息如下：");
            $section->tsv2html(
                -file   => $f_stat,
                -name   => "wgMLST 结果统计",
                -header => 0);
            `sed -i '1d' $f_stat`;
        }

        FileLink($section, "wgMLST详细结果 <link=./wgMLST.match.info.tsv>");

    }

}

sub cgsnp {
    # cgSNP 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $f_pic = "$fulldir/cgSNP.png";
    my $f_tree = "$fulldir/cgSNP.tre";

    if (-e $f_pic) {
        my $section = $report->section(id => "cgsnp");
        $section->menu("cgSNP 分析", -icon => "./src/image/cgsnp.png");
        my $desc_analysis = <<CONT;
cgSNP 分析结果如下图所示：
CONT
        $section->desc($desc_analysis);
        $section->img2html(
            -file => "cgSNP.png",
            -name => "cgSNP 分析结果"
        );
        FileLink($section, "进化树文件下载 <link=./cgSNP.tre>");
    }
}

sub wgsnp {
    # cgSNP 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $f_pic = "$fulldir/wgSNP.png";
    my $f_tree = "$fulldir/wgSNP.tre";

    if (-e $f_pic) {
        my $section = $report->section(id => "wgsnp");
        $section->menu("wgSNP 分析", -icon => "./src/image/wgsnp.png");
        my $desc_analysis = <<CONT;
wgSNP 分析结果如下图所示：
CONT
        $section->desc($desc_analysis);
        $section->img2html(
            -file => "wgSNP.png",
            -name => "wgSNP 分析结果"
        );
        FileLink($section, "进化树文件下载 <link=./wgSNP.tre>");

    }
}

sub typing {
    # 血清型分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $file = "$fulldir/typing.tsv";

    if (-e $file) {
        my $section = $report->section(id => "typing");
        $section->menu("病原微生物分型", -icon => "./src/image/typing.png");
        my $desc_analysis = <<CONT;
$blank6
CONT
        $section->desc($desc_analysis);
        $section->tsv2html(
            -file   => $file,
            -name   => "分型结果",
            -header => 1);
    }
}

sub mlva {
    # MLVA 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $file = "$fulldir/mlva.tsv";

    if (-e $file) {
        my $section = $report->section(id => "mlva");
        $section->menu("MLVA 分型分析", -icon => "./src/image/mlva.png");

        $section->tsv2html(
            -file   => $file,
            -name   => "MLVA 分型结果",
            -header => 1
        );
    }
}

sub cansnp {
    # canSNP 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $file = "$fulldir/cansnp.tsv";

    if (-e $file) {
        my $section = $report->section(id => "cansnp");
        $section->menu("canSNP 分型分析", -icon => "./src/image/cansnp.png");

        $section->tsv2html(
            -file   => $file,
            -name   => "canSNP 分型结果",
            -header => 1
        );

        FileLink($section, "canSNP详细结果 <link=./canSNP_detail.tsv>");

    }
}

sub spa {
    # Spa 分型报告
    my $report = shift;
    my $fulldir0 = $fulldir;
    my $fulldir = "$fulldir0";
    my $file = "$fulldir/Spa.tsv";

    if (-e $file) {
        my $section = $report->section(id => "spa");
        $section->menu("Spa 分型分析", -icon => "./src/image/spa.png");

        $section->tsv2html(
            -file   => $file,
            -name   => "Spa 分型结果",
            -header => 1
        );
    }
}

sub usage {
    my $help = <<HELP;

Useage:
         perl $0 -outdir Pipe

Options: -species STR     The species name [required]
         -outdir  STR     the out put dir [required]
         -nonlazy INT     lazy or non-lazy [default 1]
         -company STR     The company name [default 微远基因]
         -help|h          output this information

HELP
    print $help;
    exit 0;

}
