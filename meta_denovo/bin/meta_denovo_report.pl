my $app;
my $GDHRPATH;
BEGIN {
    use FindBin qw($Bin $RealBin);
    use YAML qw(LoadFile);
    my $f_software = "/sdbb/share/pipeline/16S/config/software.yml";
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
use Data::Dumper;
use Cwd qw(getcwd abs_path);


=head1 Description
    Genome mapping and mapped reads assembly.

=head1 Usage
    perl meta_denovo_report.pl
    -outdir <s>             The work directory for analysis, default: "./".
    -name <s>               The name of the sample id.
    -ref <s>                The name of refequence genome.
    -help                   The help message.

=head1 Exmple
    perl meta_denovo_report.pl -workdir ./ 
=cut

my ($outdir,$name,$ref,$help);
GetOptions(
    "outdir:s"=>\$outdir,
    'name:s'=>\$name,
    'ref:s'=>\$ref,
    "help"=>\$help
);
die `pod2text $0` if ($help || !$name);

my $config = "$Bin/../des_meta_denovo.txt";
my $referance_article = "$Bin/../referance_meta_denovo.txt";

$outdir ||= getcwd();
$outdir = abs_path($outdir);

my $blank6 = blank(6);

my %hash;
open A,$config || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    $hash{$a[0]} = $a[1];
}
close A;

my $ref;
open R,$referance_article || die $!;
while(<R>){
    chomp;
    $ref .= $_;
}
close R;

# 新建一个HTML对象
my $report;
$report = GDHR->new(-outdir => "$outdir",
    -pipe                   => "宏基因组数据分析报告",
    -nonlazy                => 0);

# 新建一个section对象
my $section;
$section = $report->section(id => "introduction");
$section->menu("项目概述",-icon => "src/image/gaishu.png");
# 添加描述信息
my $main = $hash{'main'};
$section->desc($main);

$section = $report->section(id => "pipeline");
$section->menu("信息分析流程",-icon => "src/image/liucheng.png");

# 添加图片
$section->img2html(
    -file => "./src/image/meta_denovo_pipe.png",
    -name => "信息分析流程图");

$section = $report->section(id => "fastq_deal");
$section->menu("数据质控",-icon => "src/image/cexuzhikong.png");


$section->desc("过滤前测序数据碱基质量分布图如下：");
my @name = ["imgs2html2"];
$section->img2html2(
    -file1 => "01.Filter/FastQC/${name}_raw_1.png",
    -name1  => "reads1",
    -file2 => "01.Filter/FastQC/${name}_raw_2.png",
    -name2  => "reads2",
    -names  => "Reads2");


$section->desc("过滤后测序数据碱基质量分布图如下：");
my @name = ["imgs2html2"];
$section->img2html2(
    -file1 => "01.Filter/FastQC/${name}_clean_1.png",
    -name1  => "reads1",
    -file2 => "01.Filter/FastQC/${name}_clean_2.png",
    -name2  => "reads2",
    -names  => "Reads2");


my $fq_qc = $hash{'fq_qc'};
$section->desc($fq_qc);

$section->tsv2html(
    -file   => "$outdir/01.Filter/${name}_qc_stat.txt",
    -name   => "过滤前后碱基信息统计表",
    );


$section = $report->section(id => "assembly");
$section->menu("基因组组装",-icon => "src/image/3.ass.png");

my $spades = $hash{'spades'};
$section->desc($spades);

$section->desc("组装结果如下表所示：");
$section->tsv2html(
    -file   => "$outdir/02.Assembly/${name}_fa.stat.txt",
    -name   => "组装结果统计表");

$section->desc("组装出的contig分布图如下图所示：");
$section->img2html(
    -file => "02.Assembly/${name}_contigs_length_distribute.png",
    -name => "Contig长度分布图");

my $fasta_dir = "02.Assembly/$name.genome.fa";
FileLink($section, "组装的基因组<link=$fasta_dir>");

$section = $report->section(id => "species_abun_taxon");
$section->menu("物种丰度及物种注释",-icon => "src/image/abun.png");


my $contig_abun_and_taxon = $hash{'contig_abun_and_taxon'};
$section->desc($contig_abun_and_taxon);

$section->tsv2html(
    -file   => "$outdir/03.Contig/${name}_Contig_abun_and_lineage.txt",
    -name   => "物种丰度和注释结果",
    );

my $abun_taxon = "03.Contig/${name}_Contig_abun_and_lineage.xlsx";
FileLink($section, "物种丰度和注释结果表完整版<link=$abun_taxon>");


$section = $report->section(id => "gene");
$section->menu("基因预测及丰度分析",-icon => "src/image/4.pre.png");
$section->submenu("基因预测");

my $gene_pre = $hash{'gene_pre'};
$section->desc($gene_pre);

my $pep_fasta_dir = "04.Gene/$name.faa";
my $gff_dir = "04.Gene/$name.gff";

FileLink($section, "基因预测gff文件<link=$gff_dir>");
FileLink($section, "基因预测蛋白序列<link=$pep_fasta_dir>");

$section->img2html(
    -file => "./04.Gene/${name}_gene_length_distribute.png",
    -name => "基因长度分布图");

$section->submenu("基因丰度分析");
my $gene_abun = $hash{'gene_abun'};
$section->desc($gene_abun);

$section->tsv2html(
    -file   => "$outdir/04.Gene/${name}_gene_abun.txt",
    -name   => "基因丰度表",
    );

my $gene_ab = "04.Gene/${name}_gene_abun.xlsx";
FileLink($section, "基因丰度表完整版<link=$gene_ab>");

$section = $report->section(id => "function");
$section->menu("功能数据库注释",-icon => "src/image/5.anno.png");
$section->submenu("GO注释");

my $go = $hash{'go'};
$section->desc($go);

my $GO_png = "$outdir/05.Function/GO/${name}_GO.png";
if (-e $GO_png){
$section->tsv2html(
    -file   => "$outdir/05.Function/GO/${name}_GO2Gene.txt",
    -name   => "GO注释表",
    -max_chars => 30,
    -header => 1);

my $go_anno = "05.Function/GO/${name}_GO2Gene.xlsx";
FileLink($section, "GO注释结果表完整版<link=$go_anno>");

$section->img2html(
    -file => "05.Function/GO/${name}_GO.png",
    -name => "GO注释");
}
else{
    $section->desc("没有基因注释到GO数据库");
}

#####################################################

$section->submenu("KEGG注释");
my $kegg = $hash{'kegg'};
$section->desc($kegg);

my $KEGG_png = "$outdir/05.Function/KEGG/${name}_KEGG.png";
if (-e $KEGG_png){
    $section->tsv2html(
        -file   => "$outdir/05.Function/KEGG/${name}_KEGG.txt",
        -name   => "KEGG注释表",
        -max_chars => 30,
        -header => 1);


    my $kegg_anno = "05.Function/KEGG/${name}_KEGG.xlsx";
    FileLink($section, "KEGG注释结果表完整版<link=$kegg_anno>");

    $section->img2html(
        -file => "05.Function/KEGG/${name}_KEGG.png",
        -name => "KEGG注释");
}

else{
    $section->desc("没有基因注释到KEGG数据库");
}

$section->submenu("COG注释");
my $cog = $hash{'cog'};
$section->desc($cog);

$section->img2html(
    -file => "05.Function/COG/${name}_COG.png",
    -name => "COG注释");

$section->submenu("CARD注释");
my $card = $hash{'card'};
$section->desc($card);

$section->tsv2html(
    -file   => "$outdir/05.Function/CARD/${name}_CARD.txt",
    -name   => "CARD注释表",);

my $card_anno = "05.Function/CARD/${name}_CARD.xlsx";
FileLink($section, "GO注释结果表完整版<link=$card_anno>");
    

$section->submenu("VFDB注释");
my $vfdb = $hash{'vfdb'};
$section->desc($vfdb);

$section->tsv2html(
    -file   => "$outdir/05.Function/VFDB/${name}_VFDB.txt",
    -name   => "VFDB注释表",);

my $vfdb_anno = "05.Function/VFDB/${name}_VFDB.xlsx";
FileLink($section, "GO注释结果表完整版<link=$vfdb_anno>");

$section = $report->section(id => "referance");
$section->menu("参考文献",-icon => "src/image/6.pub.png");

$section->desc($ref);

`cp -rfL $Bin/../report/src/image/* $outdir/src/image`;

# 添加一个
## 生成HTML报告
$report->write();

