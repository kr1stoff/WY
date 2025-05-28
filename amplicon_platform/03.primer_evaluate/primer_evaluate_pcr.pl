#!/usr/bin/perl -w 
#########################################################
# Author: stone
# Created Time : 2023.05.15 10:00:00 AM CST
# File Name: primer_evaluate_qpcr.pl
# Version: v1.0
#########################################################

use utf8;
use strict;
use warnings;
use Encode;
use FindBin qw($Bin);
use lib "$Bin/lib";
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use GACP qw(parse_config);
use Cwd qw(getcwd abs_path);
use Digest::MD5 qw(md5_hex);

=head1 Description
    This script is used to Check for specificity and inclusiveness of the primers.

=head1 Usage
    $0 -p <primer file> -o [out dir] -n [nt database] -g <genome dir> -t <name2taxid info file>
    $0 -p <primer file> -n F

=head1 Parameters
    perl meta_denovo_and_anno.pl
    -p|-primer          <s>     Input primer file with Excel format
    -b|-blast           <s>     Whether blast to nt, T|F, default:T, if -blast F, no need to set nt_database
    -n|-nt_database     <s>     The nt blast database file, default:config.txt
    -g|-genome_dir      <s>     The genome reference directory, include all.00.*, assembly_chromosome.tsv
    -t|-name2taxid      <s>     The name2taxid info file, default:config.txt
    -s|-specificity     <s>     Whether check for specificity, T|F, default:T
    -i|-inclusiveness   <s>     Whether check for inclusiveness, T|F, default:T
    -r|-run              <s>     Whether run the shell, T|F, default:T
    -e|-evalue          <s>     The blast e-value, default:800000 ,if the min primer length > 18bp, evalue can be smaller
    -o|-outdir          <s>     The outdir directory, default:'./'
    -h|-help            <s>     The help message
=cut

my ($primerfile,$nt_database,$genome_dir,,$name2taxid,$specificity,$inclusiveness,$run,$blast,$evalue,$outdir,$help);
GetOptions(
    "p|primer=s"        => \$primerfile,
    "b|blast:s"         => \$blast,
    "n|nt_database:s"   => \$nt_database,
    "g|genome_dir:s"    => \$genome_dir,
    "t|name2taxid:s"    => \$name2taxid,
    "s|specificity:s"   => \$specificity,
    "i|inclusiveness:s" => \$inclusiveness,
    "r|run:s"           => \$run,
    "e|evalue:s"        => \$evalue,
    "o|outdir:s"        =>  \$outdir,
    "h|help:s"          => \$help
    );

die `pod2text $0` if ((!$primerfile) or ($help));

# config file
my $config_file = "$Bin/config.txt";
my $nt = parse_config($config_file,"nt");
my $name2taxid_info = parse_config($config_file,"name2taxid_info");
my $taxid_dir = parse_config($config_file,"taxid_dir");
my $all_genome_dir = parse_config($config_file,"all_genome_dir");
# /home/yehui/software/python3/Python-3.8.15/bin/in2csv

#sortware
my $blastn = parse_config($config_file,"blastn");
my $taxonkit = parse_config($config_file,"taxonkit");
my $bedtools = parse_config($config_file,"bedtools");
my $in2csv = parse_config($config_file,"in2csv");

my $primer_blast_dir = "$Bin/../primer_blast";
my $conservative_dir = "$Bin/../conservative";
my $common_bin = "$Bin/../common_bin";

$outdir ||= "./";
`mkdir -p $outdir`;
$outdir = abs_path($outdir);
$primerfile = abs_path($primerfile);

$nt_database ||= $nt;
$name2taxid ||= $name2taxid_info;

$blast ||= "T";
$specificity ||= "T";
$inclusiveness ||= "T";
$run ||= "T";
$evalue ||= "800000";

my $shell_dir="$outdir/00.shell";
`mkdir $shell_dir` unless (-d $shell_dir);

my $specificity_dir = "$outdir/01.specificity";
`mkdir -p $specificity_dir`;

my $inclusiveness_dir = "$outdir/02.inclusiveness";
`mkdir -p $inclusiveness_dir`;

my $stat_dir = "$outdir/03.stat";
`mkdir -p $stat_dir`;

my @shell_list;

#格式转换，Excel文件转CSV文件
if($primerfile =~ /\.xls/){
    `csvtk xlsx2csv -t $primerfile > $outdir/Internal_Primers.txt`;
}else{
    `cat $primerfile > $outdir/Internal_Primers.txt`;
}

#创建物种名，taxid,基因组数目对应哈希
my %species2taxid;
open IN,$name2taxid or die $!;  # name to taxid info file
while(<IN>){
    chomp;
    next if /#/;
    my @line = split /\t/,$_;
    $species2taxid{$line[0]} = $line[1];
}
close IN;

open PRI,"$outdir/Internal_Primers.txt" || die $!;
open TABLE_PRIMER,">$outdir/primer.table.taxid" || die $!;
open BED,">$outdir/out.product.bed" || die $!;

my %hash_fa;
my ($taxid,$spe_name);
while(<PRI>){
    chomp;
    next if /^#/;
    my @line = split /\t/,$_;
    my ($fprimer,$rprimer,$name_tmp) = @line[1,2,3];

    # start和end都是1base
    my ($chr,$amp_start,$amp_end) = @line[4,5,6];
    my $amp_start_0base = $amp_start - 1;
    my $amp_end_1base = $amp_end + 1;

    my @name_tmps = split /-/,$name_tmp;
    $spe_name = $name_tmps[0];
    my $p_name = join ('-',@name_tmps[1..$#name_tmps]);
    if ($species2taxid{$spe_name}){
        $taxid = $species2taxid{$spe_name};
    }
    else{
        die "exp.txt文件中,物种名不对,或者物种名没有对应的taxid\n";
    }

    my $md5_fprimer = md5_hex($fprimer);
    my $md5_rprimer = md5_hex($rprimer);

    $hash_fa{$md5_fprimer} = $fprimer;
    $hash_fa{$md5_rprimer} = $rprimer;

    print TABLE_PRIMER "$p_name\t$md5_fprimer\t$md5_rprimer\t$taxid\n";
    # print BED "$chr\t$amp_start\t$amp_end_1base\n";
    print BED "$chr\t$amp_start_0base\t$amp_end\n";
}
close TABLE_PRIMER;
close BED;

print "Species: $spe_name\n";
print "Taxid: $taxid\n";

$genome_dir ||= "$all_genome_dir/$taxid";

open FA,">$outdir/primer.fa" || die $!;
foreach my $id (sort keys %hash_fa){
    print FA ">$id\n$hash_fa{$id}\n";
}
close FA;

my $outfotmat = "\'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs staxid sscinames sstrand\'";


if ($specificity eq "T"){
    my $workdir1 = "$outdir/01.specificity";
    my $cmd_spe = "echo Start specificity evaluate\n";
    $cmd_spe .= "cd $workdir1\n";
    if ($blast eq "T"){
        $cmd_spe .= "python $common_bin/blast_nt_cut.py --fasta $outdir/primer.fa --nt_database $nt_database --out $outdir/blastn.nt.out\n";
        # $cmd_spe .= "export BLASTDB=/sdbb/bioinfor/yehui/nt\n";
        # $cmd_spe .= "$blastn -task blastn-short -query $outdir/primer.fa -db $nt_database -num_threads 60 -word_size 9 -evalue $evalue -max_target_seqs 9999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt $outfotmat -out $outdir/blastn.nt.out\n";
    }
    else{
        my $st = "\n# 无需nt库比对, 直接读取 $workdir1/blastn.nt.out\n";
        print encode('UTF-8', $st);
    }
    my $taxid_list;
    if (-e "$taxid_dir/$taxid.taxid.list"){
        $taxid_list = "$taxid_dir/$taxid.taxid.list";
    }
    else{
        $cmd_spe .= "$taxonkit list -i $taxid --data-dir /home/lanlei/.taxonkit/ | awk '{print \$1}' | grep -v '^\$' > $workdir1/$taxid.taxid.list\n";
        $taxid_list = "$workdir1/$taxid.taxid.list"
    }

    $cmd_spe .= "perl $primer_blast_dir/filter_blast.pl $outdir/blastn.nt.out > $workdir1/blastn.nt.out.filter\n";
    $cmd_spe .= "mkdir -p $workdir1/01.primer\n";
    $cmd_spe .= "perl $primer_blast_dir/split_workdir.pl $outdir/primer.table.taxid $workdir1/blastn.nt.out.filter 01.primer\n";
    ###
    $cmd_spe .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $primer_blast_dir/judge_primer_blast.pl $taxid_list 01.primer/\$id/\$id.blast > 01.primer/\$id/\$id.blast.same; done\n";
    $cmd_spe .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $primer_blast_dir/same_stat.pl 01.primer/\$id/\$id.blast.same > 01.primer/\$id/\$id.blast.same.stat; done\n";
    $cmd_spe .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $primer_blast_dir/blast_primer_pair.pl 01.primer/\$id/\$id.input.table 01.primer/\$id/\$id.blast.same \$id > 01.primer/\$id/\$id.out ; done\n";
    $cmd_spe .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $primer_blast_dir/judge_primer_blast.pl $taxid_list 01.primer/\$id/\$id.out 3 > 01.primer/\$id/\$id.out.same; done\n";
    $cmd_spe .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $primer_blast_dir/same_stat.pl 01.primer/\$id/\$id.out.same > 01.primer/\$id/\$id.out.same.stat; done\n";
    $cmd_spe .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $primer_blast_dir/percent_stat.pl \$id 01.primer/\$id/\$id.input.table 01.primer/\$id/\$id.blast.same.stat 01.primer/\$id/\$id.out.same.stat > 01.primer/\$id/\$id.all.stat; done\n";
    ### 合并结果
    $cmd_spe .= "if [ -f \"primer.all.stat\" ]; then \n rm primer.all.stat \nfi\n";
    $cmd_spe .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do cat 01.primer/\$id/\$id.all.stat | grep -v '#' >> primer.all.stat; done\n";

    generateShell("$shell_dir/S1_specificity.sh", $cmd_spe ,\@shell_list);
}


if ($inclusiveness eq "T"){
    my $workdir2 = "$outdir/02.inclusiveness";
    my $cmd_inc = "echo Start inclusiveness evaluate\n";
    $cmd_inc .= "cd $workdir2\n";
    $cmd_inc .= "$blastn -task blastn-short -query $outdir/primer.fa -db $genome_dir/all -num_threads 60 -word_size 9 -evalue 1000 -max_target_seqs 99999 -qcov_hsp_perc 85 -perc_identity 85 -penalty -1 -reward 1 -outfmt $outfotmat -out $workdir2/blastn.genome.out\n";
    $cmd_inc .= "perl $primer_blast_dir/filter_blast.pl $workdir2/blastn.genome.out > $workdir2/blastn.genome.out.filter\n";
    $cmd_inc .= "mkdir -p $workdir2/01.primer\n";
    $cmd_inc .= "perl $primer_blast_dir/split_workdir.pl $outdir/primer.table.taxid $workdir2/blastn.genome.out.filter 01.primer\n";
    ###
    $cmd_inc .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $primer_blast_dir/blast_primer_pair.pl 01.primer/\$id/\$id.input.table 01.primer/\$id/\$id.blast \$id > 01.primer/\$id/\$id.out ; done\n";
    $cmd_inc .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do perl $conservative_dir/epcr_result_deal.pl $genome_dir/assembly_chromosome.tsv 01.primer/\$id/\$id.out 01.primer/\$id/\$id.assembly_chromosome.copy 01.primer/\$id/\$id.out.chr 01.primer/\$id/\$id.all.stat ; done\n";
    ### 额外的结果
    $cmd_inc .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do grep -v '#name' 01.primer/\$id/\$id.out | awk -F '\\t' '{print \$2\"\\t\"\$6-1\"\\t\"\$14\"\\t\"\$2\"\\t1\\t\"\$20}' > 01.primer/\$id/\$id.bed ; done\n";
    $cmd_inc .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do $bedtools getfasta -fi $genome_dir/all.fna -fo 01.primer/\$id/\$id.bed.fasta -bed 01.primer/\$id/\$id.bed -s; done\n";

    ### 合并结果
    $cmd_inc .= "if [ -f \"primer.all.stat\" ]; then \n rm primer.all.stat \nfi\n";
    $cmd_inc .= "cat $outdir/primer.table.taxid | while read -r id f r taxid; do cat 01.primer/\$id/\$id.all.stat | grep -v '#' >> primer.all.stat; done\n";

    generateShell("$shell_dir/S2_inclusiveness.sh", $cmd_inc ,\@shell_list);
}

if ($specificity eq "T" and $inclusiveness eq "T"){
    my $workdir3 = "$outdir/03.stat";
    my $cmd_stat = "echo Start result stat\n";
    $cmd_stat .= "cd $workdir3\n";
    $cmd_stat .= "$bedtools getfasta -fi $genome_dir/ref.fna -fo out.product.bed.fasta -bed $outdir/out.product.bed\n";
    $cmd_stat .= "perl $Bin/stat/merge_all_stat_pcr.pl $outdir/01.specificity/primer.all.stat $outdir/02.inclusiveness/primer.all.stat $outdir/Internal_Primers.txt $taxid out.product.bed.fasta\n";
    generateShell("$shell_dir/S3_stat.sh", $cmd_stat ,\@shell_list);
}

open MAIN,">$shell_dir/main.sh" || die $!;
print MAIN join (" && \\\n",@shell_list) . "\n";
close MAIN;

if ($run eq "T"){
    `bash $shell_dir/main.sh`;
}

sub generateShell{
    my ($output_shell, $content, $outshell ,$finish_string) = @_;
    unlink glob "$output_shell.*";
    $finish_string ||= "Still_waters_run_deep";
    chomp $content;
    open OUT,">$output_shell" or die "Cannot open file $output_shell:$!";
    print OUT "#!/bin/bash\n";
    print OUT "echo ==========start at : `date` ==========\n";
    print OUT "set -e \n";
    print OUT "$content && \\\n";
    print OUT "echo ==========end at : `date` ========== && \\\n";
    print OUT "echo $finish_string 1>&2 && \\\n";
    print OUT "echo $finish_string > $output_shell.sign\n";
    close OUT;
    my $run = "bash $output_shell 1>$output_shell.e 2>$output_shell.o";
    push @{$outshell}, $run;
}
