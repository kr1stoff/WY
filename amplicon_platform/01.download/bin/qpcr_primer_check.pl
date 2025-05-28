#!/usr/bin/perl
#########################################################
# Author: handsye
# Created Time : Mon 09 Jan 2023 10:41:08 AM CST
# File Name: qpcr_primer_check.pl
# Version: v0.1.0
#########################################################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils ':all';
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($INFO);

# 帮助文档
=head1 Description

    This script is used to Check for specificity and inclusiveness of the primers.

=head1 Usage

    $0 -p <primer file> -o [output dir] -n [prefix name] -b [blast database] -i <species info> -r <reference dir>

=head1 Parameters

    -p  [str]   Input primer file with Excel format
    -o  [str]   Output file to which directory(optional)
    -b  [str]   Blast database file
    -i  [str]   Input species information file
    -r  [str]   Reference file directory
    -d  [str]   Do specificity
    -h  [str]   Print this information
=cut

my ($primerfile, $blastdatabase, $allrefdir, $outdir, $info, $do, $help);
GetOptions(
    "p|primer=s" => \$primerfile,
    "b|blast=s"  => \$blastdatabase,
    "r|refdir=s" => \$allrefdir,
    "o|outdir:s" => \$outdir,
    "i|info=s"   => \$info,
    "d|do=s"     => \$do,
    "h|help:s"   => \$help
);

die `pod2text $0` if ((!$primerfile) or ($help));

$do = "T";

#创建物种名，taxid,基因组数目对应哈希
my %species2taxid;
my %species2num;
open IN, "<$info" or die $!;
while (<IN>) {
    chomp;
    next if $. == 1;
    my @line = split /\t/, $_;
    $species2taxid{$line[0]} = $line[1];
    $species2num{$line[0]} = $line[2];
}
close IN;

# #创建输出文件夹
if ($outdir) {
    # `rm -rf $outdir`;
    `mkdir -p $outdir`;
}
else {
    # `rm -rf ./primerexp-check`;
    $outdir = "./primerexp-check";
    `mkdir -p $outdir`;
}

#定义物种-引物-基因组哈希
my %taxid2primer;
#定义物种-基因组哈希表
my %taxid2genom;
#定义物种-基因组-位置哈希表
my %genom2pos;

#qPCR,内引物
INFO("Prepare");
open OUT, ">$outdir/qPCR.primer.check.xls" or die $!;
open TNGS, ">$outdir/2tNGS.bed" or die $!;
if ($do eq "T") {
    print OUT "ID\tLeft_Seq\tRight_Seq\tProbe_Seq\tAmplicon_Name\tChrom\tAmp_Start\tAmp_End\tAmp_Length\tLeft_Tm\tRight_Tm\tProbe_Tm\tLeft_Length\tRight_Length\tProbe_Length\tLeft_GC\tRight_GC\tProbe_GC\tTaxid\tAmplicon_Seq\t1-specificity\t2-specificity\tprimer_inclusiveness\tprobe-specificity\tprobe_inclusiveness\tuse\n";
}
else {
    print OUT "ID\tLeft_Seq\tRight_Seq\tProbe_Seq\tAmplicon_Name\tChrom\tAmp_Start\tAmp_End\tAmp_Length\tLeft_Tm\tRight_Tm\tProbe_Tm\tLeft_Length\tRight_Length\tProbe_Length\tLeft_GC\tRight_GC\tProbe_GC\tTaxid\tAmplicon_Seq\tprimer_inclusiveness\tprobe_inclusiveness\tuse\n";
}
`dos2unix $primerfile`;

my $refseq_num; #GCF基因组数量
my $all_num;    #GCF+GCA基因组数量
my $reftaxid;
#比对NT数据库,比对物种基因组数据库
open NT, ">align2nt.fa" or die $!;
open WZ, ">align2wz.fa" or die $!;
open PRI, "<$primerfile" or die $!;
while (<PRI>) {
    chomp;
    next if $_ =~ /^#/;
    my @line = split /\t/, $_;
    my @line2 = split /-/, $line[27];
    $line[27] =~ s/ /-/g;
    $line[27] =~ s/\[//g;
    $line[27] =~ s/\]//g;

    #引物序列1
    my $F = join("-", $line[27], "F");
    print NT ">$F\n$line[6]\n";

    #引物序列2
    my $R = join("-", $line[27], "R");
    print NT ">$R\n$line[12]\n";

    #探针序列
    my $P = join("-", $line[27], "P");
    print NT ">$P\n$line[18]\n";
    print WZ ">$P\n$line[18]\n";

    if ($. == 2) {
        my $species = $line2[0];
        if ($species2taxid{$species}) {
            $reftaxid = $species2taxid{$species};
        }
        if (-e "$allrefdir/$reftaxid/assembly_chromosome.tsv") {
            `grep "GCF" $allrefdir/$reftaxid/assembly_chromosome.tsv > tmp.assembly_chromosome.tsv`;
            $refseq_num = `wc -l tmp.assembly_chromosome.tsv|cut -f1 -d " "`;
            $all_num = `wc -l $allrefdir/$reftaxid/assembly_chromosome.tsv|cut -f1 -d " "`;
            $all_num = $all_num - 1;
        }
        else {
            print "assembly_chromosome.tsv数据表不存在，请检查！\n";
        }
    }

}
close NT;
close PRI;
close WZ;

# RCB: 暂时剔除
INFO("Blast against the nt database");
# &blast("align2nt.fa", $blastdatabase);

# if ($refseq_num >= 10) {
#     if (-e "$allrefdir/$reftaxid/refseq.fna") {
#         &blast2("align2wz.fa", "$allrefdir/$reftaxid/refseq");
#     }
#     else {
#         print "refseq.fna数据库不存在，请检查！\n";
#     }
# }
# else {
#     if (-e "$allrefdir/$reftaxid/all.fna") {
#         &blast2("align2wz.fa", "$allrefdir/$reftaxid/all");
#     }
#     else {
#         print "all.fna数据库不存在，请检查！\n";
#     }
# }

INFO("Parse the blast result");
open PRI, "<$primerfile" or die $!;
open NORUN, ">$outdir/norun-qPCR-Primers.csv" or die $!; #输出有问题的引物对
while (<PRI>) {
    chomp;
    next if $_ =~ /^#/;
    my @line = split /\t/, $_;
    my $amp_start = $line[25] + $line[1] - 1;
    my $amp_end = $line[25] + $line[7] - 1;
    my $y; #primer express 设计结果
    $y = join("\t", $line[6], $line[12], $line[18], $line[27], $line[24], $amp_start, $amp_end, $line[22], $line[4], $line[10], $line[16], $line[3], $line[9], $line[15], $line[5], $line[11], $line[17]);
    my $primer_probe = join("-", $line[6], $line[12], $line[18]);
    my $primer = join("-", $line[6], $line[12]);
    my @line2 = split /-/, $line[27];
    my $species = $line2[0];
    my $reftaxid;     #物种taxid
    my $id;           #引物唯一ID
    my $amplicon_seq; #扩增子序列

    my ($tp1, $tp2, $tp3, $bao1, $bao2);

    if ($species2taxid{$species}) {
        $reftaxid = $species2taxid{$species};
        DEBUG
        ("$species\t$reftaxid\n");
        `echo $primer_probe $reftaxid > id.txt`;
        $id = &md5("id.txt");
        `rm id.txt`;

        #创建物种包含的所有子taxid的对应表
        # `taxonkit list --show-rank --show-name --ids $reftaxid --data-dir /home/yehui/software/TaxonKit > idtemp.txt`;
        my %taxid2species;
        open IN, "<idtemp.txt" or die $!;
        while (<IN>) {
            chomp;
            next if $_ =~ /^$/; #跳过空行
            $_ =~ s/^\s*//;
            my @line = split / /, $_;
            $taxid2species{$line[0]} = $species;
        }
        close IN;
        #print Dumper(\%taxid2species);

        #创建taxid，基因组名，染色体名对应表
        my %chrom2genom;
        my $gcf_num;     #GCF基因组数量
        my $gcf_gca_num; #GCF+GCA基因组数量
        if (-e "$allrefdir/$reftaxid/assembly_chromosome.tsv") {
            my $file;
            `grep "GCF" $allrefdir/$reftaxid/assembly_chromosome.tsv > tmp.assembly_chromosome.tsv`;
            $gcf_num = `wc -l tmp.assembly_chromosome.tsv|cut -f1 -d " "`;
            $gcf_gca_num = `wc -l $allrefdir/$reftaxid/assembly_chromosome.tsv|cut -f1 -d " "`;
            $gcf_gca_num = $gcf_gca_num - 1;
            if ($gcf_num < 10) {
                `sed '1d' "$allrefdir/$reftaxid/assembly_chromosome.tsv" > tmp2.assembly_chromosome.tsv`;
                $file = "tmp2.assembly_chromosome.tsv";
            }
            else {
                $file = "tmp.assembly_chromosome.tsv";
            }
            #open IN,"<$allrefdir/$reftaxid/tmp.assembly_chromosome.tsv" or die $!;
            #open IN,"<tmp.assembly_chromosome.tsv" or die $!;
            open IN, "<$file" or die $!;
            while (<IN>) {
                chomp;
                #next if $.==1;
                my @line = split /\t/, $_;
                my @line2 = split /,/, $line[2];
                foreach my $i (@line2) {
                    $chrom2genom{$i} = $line[1];
                }
                if (grep {$_ eq $line[1]} @{$taxid2genom{$reftaxid}}) {
                    next;
                }
                else {
                    push @{$taxid2genom{$reftaxid}}, $line[1];
                }
                $genom2pos{$reftaxid}{$line[1]} = $.;
            }
            close IN;
        }
        else {
            print "assembly_chromosome.tsv数据表不存在，请检查！\n";
            %chrom2genom = ();
        }
        #`rm tmp.assembly_chromosome.tsv`;


        #引物序列1特异性
        my $namename = $line[27];
        $namename =~ s/ /-/g;
        $namename =~ s/\[//g;
        $namename =~ s/\]//g;
        my $F = join("-", $namename, "F");
        `grep -w $F out.blast > primer_F.blast`;
        $tp1 = &parseblast("primer_F.blast", $species, \%taxid2species);

        #引物序列2特异性
        my $R = join("-", $namename, "R");
        `grep -w $R out.blast > primer_R.blast`;
        $tp2 = &parseblast("primer_R.blast", $species, \%taxid2species);


        #扩增子序列（包含引物）
        open OU3, ">amp.bed" or die $!;
        print OU3 "$line[24]\t$amp_start\t$amp_end\n";
        close OU3;
        &getfa("$allrefdir/$reftaxid/ref.fna", "amp.bed");
        $amplicon_seq = `sed -n '2p' out.fa`;
        $amplicon_seq =~ s/\n//;

        #引物包容性分析
        if ($gcf_num >= 10) {
            if (-e "$allrefdir/$reftaxid/refseq.fna" && %chrom2genom) {
                &getamplicon("$allrefdir/$reftaxid/refseq.fna", $line[6], $line[12]);
                my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                #$bao1 = $len / $species2num{$species} * 100;
                $bao1 = $len / $gcf_num * 100;
                $bao1 = sprintf "%.3f", $bao1;

                my @unigenom = @$uni;
                foreach my $i (@unigenom) {
                    push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                }
            }
            else {
                print "refseq.fna数据库不存在，请检查！\n";
                $bao1 = 0;
            }
        }
        else {
            if (-e "$allrefdir/$reftaxid/all.fna" && %chrom2genom) {
                &getamplicon("$allrefdir/$reftaxid/all.fna", $line[6], $line[12]);
                my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                $bao1 = $len / $gcf_gca_num * 100;
                $bao1 = sprintf "%.3f", $bao1;

                my @unigenom = @$uni;
                foreach my $i (@unigenom) {
                    push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                }
            }
            else {
                print "all.fna数据库不存在，请检查！\n";
                $bao1 = 0;
            }
        }

        #探针特异性分析
        my $P = join("-", $namename, "P");
        `grep -w $P out.blast > primer_P.blast`;
        $tp3 = &parseblast("primer_P.blast", $species, \%taxid2species);

        #探针包容性分析
        `grep -w $P probe.blast > probe_P.blast`;
        if ($gcf_num >= 10) {
            $bao2 = &parseblast2("probe_P.blast", $gcf_num, \%chrom2genom);

        }
        else {
            $bao2 = &parseblast2("probe_P.blast", $gcf_gca_num, \%chrom2genom);
        }

    }
    else {
        ERROR("物种拉丁名错误，请检查！\n");
        $reftaxid = "xxx";
        DEBUG
        ("$species\t$reftaxid\n");
        print NORUN "$_\n";

        $tp1 = 0;
        $tp2 = 0;
        $tp3 = 0;
        $bao1 = 0;
        $bao2 = 0;
    }

    #判断引物对是否可用 
    my $use;
    if (($tp1 >= 80 || $tp2 >= 80) && $bao1 >= 80 && $bao2 >= 80) {
        $use = "YES";
        print TNGS "$line[24]\t$amp_start\t$amp_end\t$line[27]\t1\tB\tW\n";
    }
    else {
        $use = "NO";
    }

    # if($use eq "YES"){
    #     print TNGS "$line[24]\t$amp_start\t$amp_end\t$line[27]\t1\tB\tW\n";
    # }else{
    #     next;
    # }

    #输出结果
    if ($do eq "T") {
        DEBUG
        ("$line[6]\t$line[12]\t$tp1\t$tp2\t$bao1\t$tp3\t$bao2\t$use\n");
        print OUT "$id\t$y\t$reftaxid\t$amplicon_seq\t$tp1\t$tp2\t$bao1\t$tp3\t$bao2\t$use\n";
    }
    else {
        DEBUG
        ("$line[6]\t$line[12]\t$bao1\t$bao2\t$use\n");
        print OUT "$id\t$y\t$reftaxid\t$amplicon_seq\t$bao1\t$bao2\t$use\n";
    }
}
close OUT;
close PRI;
close NORUN;
close TNGS;


#输出包容性表格
INFO("Output the result");
my $landir = "$outdir/2lan";
`mkdir -p $landir`;
foreach my $x (keys %taxid2primer) {
    #taxid
    open TAX, ">$landir/$x.xls" or die $!;
    my $num = scalar(@{$taxid2genom{$x}});
    print "物种基因组数量：$num\n";
    my @zero = (0) x $num;
    my $a = join("\t", @{$taxid2genom{$x}});
    print TAX "primer\t$a\n";
    foreach my $y (keys %{$taxid2primer{$x}}) {
        #引物
        foreach my $z (@{$taxid2primer{$x}{$y}}) { #基因组
            if (exists $genom2pos{$x}{$z}) {
                $zero[$genom2pos{$x}{$z}] = 1;
            }
        }
        my $b = join("\t", @zero);
        print TAX "$y\t$b\n";

    }
    close TAX;
}


#=======================子函数==========================#
sub blast {
    #$_[0]：需要blast的fa序列
    #$_[1]：blast用到的数据库
    #`blastn -task blastn-short -query $_[0] -db $_[1] -num_threads 50 -word_size 7 -outfmt "7 staxid ssciname qseqid sseqid length qlen slen sstart send qstart qend pident nident evalue bitscore" -out out.blast`;
    `blastn -task blastn-short -query $_[0] -db $_[1] -num_threads 100 -word_size 7 -max_target_seqs 9999 -outfmt "6 staxid ssciname qlen mismatch qstart qend pident nident evalue bitscore qseqid" -out out.blast`;
}

sub blast2 {
    #$_[0]：需要blast的fa序列
    #$_[1]：blast用到的数据库
    #`blastn -task blastn-short -query $_[0] -db $_[1] -num_threads 50 -word_size 7 -outfmt "7 staxid ssciname qseqid sseqid length qlen slen sstart send qstart qend pident nident evalue bitscore" -out out.blast`;
    #`blastn -task blastn-short -query $_[0] -db $_[1] -num_threads 100 -word_size 7 -outfmt "6 staxid ssciname qlen mismatch qstart qend pident nident evalue bitscore" -out out.blast`;
    `blastn -task blastn-short -query $_[0] -out probe.blast -num_threads 100 -word_size 7 -max_target_seqs 9999 -db $_[1] -outfmt '6 qaccver saccver qlen mismatch qstart qend pident nident evalue bitscore qseqid' `;
}

sub parseblast {
    #$_[0]：blast得到的结果文件
    #$_[1]：物种名
    #$_[2]：物种包含的所有的taxid的信息，哈希表
    my %num;       #每个物种blast的基因组数量
    my %tax;       #每个物种的taxid
    my $hits = 0;  #blast结果总基因组数
    my $count = 0; #blast结果属于本物种的总数量
    my $list3 = $_[2];
    my %taxid2species = %$list3;
    #print Dumper(\%taxid2species);

    open IN, "<$_[0]" or die $!;
    while (<IN>) {
        chomp;
        my @line = split /\t/, $_;
        #过滤blast结果，根据序列最后一个碱基必须比对上、identify、mismatch三个条件
        if (($line[2] == $line[5]) && $line[6] >= 95 && $line[3] <= 2) {
            $hits += 1;
            if (exists $taxid2species{$line[0]}) {
                $count += 1;
            }
            else {
                next;
            }
        }
        else {
            next;
        }
    }
    close IN;

    #定植菌判断
    #留空
    #留空

    #print "$hits\n";
    my $p;
    if ($hits gt 0) {
        #print "$num{$_[1]}\n";
        #$p = $num{$_[1]} / $hits * 100;#特异性占比
        $p = $count / $hits * 100; #特异性占比
    }
    else {
        DEBUG
        ("blast无结果\t");
        DEBUG
        ("$hits\n");
        $p = "0";
    }

    $p = sprintf "%.3f", $p;
    DEBUG
    ("$p\n");
    return ($p);
}

sub parseblast2 {
    #$_[0]：blast得到的结果文件
    #$_[1]：物种包含的基因组数目
    #$_[2]：染色体名和基因组名对应信息，哈希表
    my $list3 = $_[2];
    my %chrom2genom = %$list3;
    my @genom; #blast结果总基因组
    open IN, "<$_[0]" or die $!;
    while (<IN>) {
        chomp;
        my @line = split /\t/, $_;
        #过滤blast结果，根据序列最后一个碱基必须比对上、identify、mismatch三个条件
        if (($line[2] == $line[5]) && $line[6] >= 95 && $line[3] <= 1) {
            if (exists $chrom2genom{$line[1]}) {
                push @genom, $chrom2genom{$line[1]};
            }
            else {
                next;
            }
        }
        else {
            next;
        }
    }
    close IN;

    #定植菌判断
    #留空
    #留空

    #print "$hits\n";
    my $p;
    my $count = $_[1];
    DEBUG
    ("$count\n");
    my @uni = uniq(@genom);

    my $len = @uni;
    DEBUG
    ("$len\n");
    if ($len gt 0) {
        #print "$num{$_[1]}\n";
        $p = $len / $count * 100; #包容性占比
    }
    else {
        print "blast无结果或过滤条件严格，请检查！！！\t";
        $p = "0";
    }

    $p = sprintf "%.3f", $p;
    DEBUG
    ("$p\n");
    return ($p);
}

sub getfa {
    #$_[0]：参考基因组
    #$_[1]：bed区间文件
    `/home/chenwenjing/pipeline/IDseq_v2.0/bin/bedtools getfasta -fi $_[0] -bed $_[1] -fo out.fa`;
}

sub getamplicon {
    #$_[0]：全部基因组
    #$_[1]：F引物
    #$_[2]：R引物
    `/home/yehui/software/bin/seqkit amplicon $_[0] -F $_[1]  -R $_[2]  -m 2 -j 20 --bed --quiet >amplicon.txt`; #设置错误匹配2
}

sub parseamplicon {
    #$_[0]：扩增片段文件
    #$_[1]：染色体名和基因组名对应信息，哈希表
    my $list3 = $_[1];
    my %chrom2genom = %$list3;
    my @genom;
    open IN, "<$_[0]" or die $!;
    while (<IN>) {
        chomp;
        my @line = split /\t/, $_;
        if (exists $chrom2genom{$line[0]}) {
            push @genom, $chrom2genom{$line[0]};
        }
        else {
            next;
        }
    }
    close IN;
    my @uni = uniq(@genom);
    my $len = @uni;
    return ($len, \@uni);
}

sub getamplicon2 {
    #$_[0]：引物fa文件
    #$_[1]：全部基因组
    `mfeprimer index $_[1]`;
    `mfeprimer -i $_[0] -d $_[1] -c 20 --fasta -o mfe`; #3'端不能有错误匹配
}

sub md5 {
    #$_[0]：引物+物种信息文件
    `md5sum $_[0] > id.md5`;
    my $id;
    open IN, "<id.md5" or die $!;
    while (<IN>) {
        chomp;
        my @line = split /\s/, $_;
        $id = $line[0];
    }
    close IN;
    `rm id.md5`;
    return ($id);
}

sub align_host{

}