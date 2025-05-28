#!/usr/bin/perl
#########################################################
# Author: handsye
# Created Time : Wed 26 Oct 2022 10:41:42 AM CST
# File Name: pcr_primer_check.pl
# Version: v0.2.0
#########################################################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils ':all';

# 帮助文档
=head1 Description

    This script is used to Check for specificity and inclusiveness of the primers.

=head1 Usage

    $0 -p <primer file> -o [output dir] -n [prefix name] -b <blast database> -i <species info> -r <reference dir> -t <qPCR|tNGS|Both>

=head1 Parameters

    -p  [str]   Input primer file with Excel format
    -o  [str]   Output file to which directory(optional)
    -n  [str]   Prefix name of the result file(optional)
    -b  [str]   Blast database file
    -i  [str]   Input species information file
    -r  [str]   Reference file directory
    -t  [str]   Type of primers
    -h  [str]   Print this information
=cut

my ($primerfile, $blastdatabase, $allrefdir, $outdir, $runname, $info, $type, $help);
GetOptions(
    "p|primer=s" => \$primerfile,
    "b|blast=s"  => \$blastdatabase,
    "r|refdir=s" => \$allrefdir,
    "o|outdir:s" => \$outdir,
    "n|name:s"   => \$runname,
    "i|info=s"   => \$info,
    "t|type=s"   => \$type,
    "h|help:s"   => \$help
);

die `pod2text $0` if ((!$primerfile) or ($help));

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

#创建输出文件夹
if ($outdir) {
    `rm -rf $outdir`;
    `mkdir -p $outdir`;
}
else {
    `rm -rf ./primergo-check`;
    $outdir = "./primergo-check";
    `mkdir -p $outdir`;
}

#格式转换，Excel文件转CSV文件
#自定义文件转特定文件名
if ($primerfile =~ /\.xls/) {
    `/home/yehui/software/python3/Python-3.8.15/bin/in2csv --sheet WY_primergo_Internal_Primers $primerfile >$outdir/Internal_Primers.csv`;
    `/home/yehui/software/python3/Python-3.8.15/bin/in2csv --sheet WY_primergo_External_Primers $primerfile >$outdir/External_Primers.csv`;
}
else {
    if (-e "$outdir/Internal_Primers.csv") {
        `mv $outdir/Internal_Primers.csv $outdir/Internal_Primers.csv.bak`;
        `cat $primerfile > $outdir/Internal_Primers.csv`;
    }
    else {
        `cat $primerfile > $outdir/Internal_Primers.csv`;
    }

    if (-e "$outdir/External_Primers.csv") {
        `mv $outdir/External_Primers.csv $outdir/External_Primers.csv.bak`;
        `cat $primerfile > $outdir/External_Primers.csv`;
    }
    else {
        `cat $primerfile > $outdir/External_Primers.csv`;
    }
}

#定义物种-引物-基因组哈希
my %taxid2primer;
#定义物种-基因组哈希表
my %taxid2genom;
#定义物种-基因组-位置哈希表
my %genom2pos;

##主流程，特异性分析和包容性分析
#tNGS,外引物
if ($type eq "tNGS") {
    if ($runname) {
        open OUT, ">$outdir/$runname.tNGS.primer.check.xls" or die $!;
    }
    else {
        open OUT, ">$outdir/tNGS.primer.check.xls" or die $!;
    }
    open PRI, "<$outdir/Internal_Primers.csv" or die $!;
    open NORUN, ">$outdir/norun-tNGS-Primers.csv" or die $!; #输出有问题的引物对
    while (<PRI>) {
        chomp;
        #next if $.==1;
        if ($. == 1) {
            my @line = split /,/, $_;
            $line[0] = "ID";
            my $x = join("\t", @line);
            print OUT "$x\t";
            print OUT "Taxid\tAmplicon_Seq\t1-specificity\t2-specificity\tSE751-specificity\tSE752-specificity\tSE1001-specificity\tSE1002-specificity\tinclusiveness\tuse\n";
            print NORUN "$_\n";
        }
        next if $. == 1;
        #print "$_\n";
        my @line = split /,/, $_;
        for (my $a = 0; $a <= $#line; $a++) {
            $line[$a] =~ s/\.0//;
        }
        my $y; #引物设计数据
        my $length = (scalar @line) - 1;
        #print "$length\n";
        if ($length == 17) {
            $y = join("\t", @line[1 .. 17]);
        }
        else {
            $y = join("\t", @line[1 .. 16], "-");
        }
        #my $y = join("\t",@line);
        my $primer = join("-", @line[1 .. 2]); #引物对信息
        # my @line2 = split/_/,$line[3];
        # my $species = join(" ",@line2[0..$#line2-2]);
        my @line2 = split /-/, $line[3];
        my $species = $line2[0];
        my $reftaxid;     #物种taxid
        my $id;           #引物唯一ID
        my $amplicon_seq; #扩增子序列
        my ($tp1, $tp2, $tp3, $tp4, $tp5, $tp6, $bao);

        if ($species2taxid{$species}) {
            $reftaxid = $species2taxid{$species};
            print "$species\t$reftaxid\n";
            `echo  $primer $reftaxid > id.txt`;
            $id = &md5("id.txt");
            `rm id.txt`;

            #创建物种包含的所有子taxid的对应表
            `taxonkit list --show-rank --show-name --ids $reftaxid --data-dir /home/yehui/software/TaxonKit > idtemp.txt`;
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
            `rm idtemp.txt`;
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
                if ($gcf_num <= 2) {
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

            #引物序列1
            open OU1, ">primer1.fa" or die $!;
            print OU1 ">$species\n$line[1]";
            close OU1;
            &blast("primer1.fa", $blastdatabase);
            $tp1 = &parseblast("out.blast", $species, \%taxid2species);

            #引物序列2
            open OU2, ">primer2.fa" or die $!;
            print OU2 ">$species\n$line[2]";
            close OU2;
            &blast("primer2.fa", $blastdatabase);
            $tp2 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序75bp一端
            open OU3, ">SE751.bed" or die $!;
            my $add1 = 75 - $line[12];
            my $SE751_end = $line[5] + $add1;
            print OU3 "$line[4]\t$line[5]\t$SE751_end\n";
            close OU3;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE751.bed");
            &blast("out.fa", $blastdatabase);
            $tp3 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序75bp另一端
            open OU4, ">SE752.bed" or die $!;
            my $add2 = 75 - $line[13];
            my $SE752_start = $line[6] - $add2;
            print OU4 "$line[4]\t$SE752_start\t$line[6]\n";
            close OU4;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE752.bed");
            &blast("out.fa", $blastdatabase);
            $tp4 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序100bp一端
            open OU5, ">SE1001.bed" or die $!;
            my $add3 = 100 - $line[12];
            my $SE1001_end = $line[5] + $add3;
            print OU5 "$line[4]\t$line[5]\t$SE1001_end\n";
            close OU5;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE1001.bed");
            &blast("out.fa", $blastdatabase);
            $tp5 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序100bp另一端
            open OU6, ">SE1002.bed" or die $!;
            my $add4 = 100 - $line[13];
            my $SE1002_start = $line[6] - $add4;
            print OU6 "$line[4]\t$SE1002_start\t$line[6]\n";
            close OU6;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE1002.bed");
            &blast("out.fa", $blastdatabase);
            $tp6 = &parseblast("out.blast", $species, \%taxid2species);

            #扩增子序列（包含引物）
            open OU7, ">amp.bed" or die $!;
            print OU7 "$line[4]\t$line[5]\t$line[6]\n";
            close OU7;
            &getfa("$allrefdir/$reftaxid/ref.fna", "amp.bed");
            $amplicon_seq = `sed -n '2p' out.fa`;
            $amplicon_seq =~ s/\n//;

            #留空
            #留空

            #引物包容性分析
            if ($gcf_num >= 10) {
                if (-e "$allrefdir/$reftaxid/refseq.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/refseq.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $species2num{$species} * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "refseq.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
            else {
                if (-e "$allrefdir/$reftaxid/all.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/all.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $gcf_gca_num * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "all.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
        }
        else {
            print "物种拉丁名错误，请检查！\n";
            $reftaxid = "xxx";
            print "$species\t$reftaxid\n";
            print NORUN "$_\n";

            $tp1 = 0;
            $tp2 = 0;
            $tp3 = 0;
            $tp4 = 0;
            $tp5 = 0;
            $tp6 = 0;
            $bao = 0;
        }

        #判断引物对是否可用
        my $use;
        if (($tp1 >= 80 || $tp2 >= 80 || $tp3 >= 80 || $tp4 >= 80 || $tp5 >= 80 || $tp6 >= 80) && $bao >= 80) {
            $use = "YES";
        }
        else {
            $use = "NO";
        }

        #输出结果
        print "$line[1]\t$line[2]\t$tp1\t$tp2\t$tp3\t$tp4\t$tp5\t$tp6\t$bao\t$use\n";
        print OUT "$id\t$y\t$reftaxid\t$amplicon_seq\t$tp1\t$tp2\t$tp3\t$tp4\t$tp5\t$tp6\t$bao\t$use\n";
    }
    close OUT;
    close PRI;
    close NORUN;

    #输出包容性表格
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
}

#巢式PCR外引物
if ($type eq "caoshi") {
    `rm -rf Blast_result`;
    `mkdir Blast_result`;
    if ($runname) {
        open OUT, ">$outdir/$runname.tNGS.primer.check.xls" or die $!;
    }
    else {
        open OUT, ">$outdir/tNGS.primer.check.xls" or die $!;
    }
    open PRI, "<$outdir/Internal_Primers.csv" or die $!;
    open NORUN, ">$outdir/norun-tNGS-Primers.csv" or die $!; #输出有问题的引物对
    while (<PRI>) {
        chomp;
        #next if $.==1;
        if ($. == 1) {
            my @line = split /,/, $_;
            $line[0] = "ID";
            my $x = join("\t", @line);
            print OUT "$x\t";
            print OUT "Taxid\tAmplicon_Seq\t1-specificity\t2-specificity\tinclusiveness\tuse\n";
            print NORUN "$_\n";
        }
        next if $. == 1;
        #print "$_\n";
        my @line = split /,/, $_;
        for (my $a = 0; $a <= $#line; $a++) {
            $line[$a] =~ s/\.0//;
        }
        my $y; #引物设计数据
        my $length = (scalar @line) - 1;
        #print "$length\n";
        if ($length == 17) {
            $y = join("\t", @line[1 .. 17]);
        }
        else {
            $y = join("\t", @line[1 .. 16], "-");
        }
        #my $y = join("\t",@line);
        my $primer = join("-", @line[1 .. 2]); #引物对信息
        # my @line2 = split/_/,$line[3];
        # my $species = join(" ",@line2[0..$#line2-2]);
        my @line2 = split /-/, $line[3];
        my $species = $line2[0];
        my $tmpid = join("-", @line2[1 .. $#line2]);
        my $reftaxid;     #物种taxid
        my $id;           #引物唯一ID
        my $amplicon_seq; #扩增子序列
        my ($tp1, $tp2, $tp3, $tp4, $tp5, $tp6, $bao);

        if ($species2taxid{$species}) {
            $reftaxid = $species2taxid{$species};
            print "$species\t$reftaxid\n";
            `echo  $primer $reftaxid > id.txt`;
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
            # `rm idtemp.txt`;
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

            #引物序列1
            open OU1, ">F.fa" or die $!;
            print OU1 ">$species\n$line[1]";
            close OU1;
            &blast("F.fa", $blastdatabase);
            `mkdir -p Blast_result/$tmpid`;
            $tp1 = &parseblast("out.blast", $species, \%taxid2species);
            `mv F.fa Blast_result/$tmpid`;
            `mv out.blast Blast_result/$tmpid/F.blast`;

            #引物序列2
            open OU2, ">R.fa" or die $!;
            print OU2 ">$species\n$line[2]";
            close OU2;
            &blast("R.fa", $blastdatabase);
            $tp2 = &parseblast("out.blast", $species, \%taxid2species);
            `mv R.fa Blast_result/$tmpid`;
            `mv out.blast Blast_result/$tmpid/R.blast`;

            #扩增子序列（包含引物）
            open OU7, ">amp.bed" or die $!;
            print OU7 "$line[4]\t$line[5]\t$line[6]\n";
            close OU7;
            &getfa("$allrefdir/$reftaxid/ref.fna", "amp.bed");
            $amplicon_seq = `sed -n '2p' out.fa`;
            $amplicon_seq =~ s/\n//;
            `mv amp.bed Blast_result/$tmpid`;

            #留空
            #留空

            #引物包容性分析
            if ($gcf_num >= 10) {
                if (-e "$allrefdir/$reftaxid/refseq.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/refseq.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    `mv amplicon.txt Blast_result/$tmpid`;
                    $bao = $len / $gcf_num * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "refseq.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
            else {
                if (-e "$allrefdir/$reftaxid/all.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/all.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    `mv amplicon.txt Blast_result/$tmpid`;
                    $bao = $len / $gcf_gca_num * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "all.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
        }
        else {
            print "物种拉丁名错误，请检查！\n";
            $reftaxid = "xxx";
            print "$species\t$reftaxid\n";
            print NORUN "$_\n";

            $tp1 = 0;
            $tp2 = 0;
            $bao = 0;
        }

        #判断引物对是否可用
        my $use;
        if (($tp1 >= 80 || $tp2 >= 80) && $bao >= 80) {
            $use = "YES";
        }
        else {
            $use = "NO";
        }

        #输出结果
        print "$line[1]\t$line[2]\t$tp1\t$tp2\t$bao\t$use\n";
        print OUT "$id\t$y\t$reftaxid\t$amplicon_seq\t$tp1\t$tp2\t$bao\t$use\n";
    }
    close OUT;
    close PRI;
    close NORUN;

    #输出包容性表格
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
}

#qPCR,内引物
if ($type eq "qPCR") {
    if ($runname) {
        open OUT, ">$outdir/$runname.qPCR.primer.check.xls" or die $!;
    }
    else {
        open OUT, ">$outdir/qPCR.primer.check.xls" or die $!;
    }
    open PRI, "<$outdir/Internal_Primers.csv" or die $!;
    open NORUN, ">$outdir/norun-qPCR-Primers.csv" or die $!; #输出有问题的引物对
    while (<PRI>) {
        chomp;
        #next if $.==1;
        if ($. == 1) {
            my @line = split /,/, $_;
            $line[0] = "ID";
            my $x = join("\t", @line);
            print OUT "$x\t";
            print OUT "Taxid\tAmplicon_Seq\t1-specificity\t2-specificity\tamplicon\tinclusiveness\tuse\n";
            print NORUN "$_\n";
        }
        next if $. == 1;
        #print "$_\n";
        my @line = split /,/, $_;
        for (my $a = 0; $a <= $#line; $a++) {
            $line[$a] =~ s/\.0//;
        }
        my $y; #引物设计数据
        my $length = (scalar @line) - 1;
        if ($length == 17) {
            $y = join("\t", @line[1 .. 17]);
        }
        else {
            $y = join("\t", @line[1 .. 16], "-");
        }
        #my $y = join("\t",@line);
        my $primer = join("-", @line[1 .. 2]);
        # my @line2 = split/_/,$line[3];
        # my $species = join(" ",@line2[0..$#line2-2]);
        my @line2 = split /-/, $line[3];
        my $species = $line2[0];
        my $reftaxid;     #物种taxid
        my $id;           #引物唯一ID
        my $amplicon_seq; #扩增子序列

        my ($tp1, $tp2, $tp3, $bao);

        if ($species2taxid{$species}) {
            $reftaxid = $species2taxid{$species};
            print "$species\t$reftaxid\n";
            `echo $primer $reftaxid > id.txt`;
            $id = &md5("id.txt");
            `rm id.txt`;

            #创建物种包含的所有子taxid的对应表
            `taxonkit list --show-rank --show-name --ids $reftaxid --data-dir /home/yehui/software/TaxonKit > idtemp.txt`;
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
                if ($gcf_num <= 2) {
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

            #引物序列1
            open OU1, ">primer1.fa" or die $!;
            print OU1 ">$species\n$line[1]";
            close OU1;
            &blast("primer1.fa", $blastdatabase);
            $tp1 = &parseblast("out.blast", $species, \%taxid2species);

            #引物序列2
            open OU2, ">primer2.fa" or die $!;
            print OU2 ">$species\n$line[2]";
            close OU2;
            &blast("primer2.fa", $blastdatabase);
            $tp2 = &parseblast("out.blast", $species, \%taxid2species);

            #插入片段
            open OU3, ">amplicon.bed" or die $!;
            print OU3 "$line[4]\t$line[8]\t$line[9]\n";
            close OU3;
            &getfa("$allrefdir/$reftaxid/ref.fna", "amplicon.bed");
            &blast("out.fa", $blastdatabase);
            $tp3 = &parseblast("out.blast", $species, \%taxid2species);

            #扩增子序列（包含引物）
            open OU4, ">amp.bed" or die $!;
            print OU4 "$line[4]\t$line[5]\t$line[6]\n";
            close OU4;
            &getfa("$allrefdir/$reftaxid/ref.fna", "amp.bed");
            $amplicon_seq = `sed -n '2p' out.fa`;
            $amplicon_seq =~ s/\n//;

            #引物包容性分析
            if ($gcf_num >= 10) {
                if (-e "$allrefdir/$reftaxid/refseq.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/refseq.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $species2num{$species} * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "refseq.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
            else {
                if (-e "$allrefdir/$reftaxid/all.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/all.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $gcf_gca_num * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "all.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
        }
        else {
            print "物种拉丁名错误，请检查！\n";
            $reftaxid = "xxx";
            print "$species\t$reftaxid\n";
            print NORUN "$_\n";

            $tp1 = 0;
            $tp2 = 0;
            $tp3 = 0;
            $bao = 0;
        }

        #判断引物对是否可用
        my $use;
        if (($tp1 >= 80 || $tp2 >= 80 || $tp3 >= 80) && $bao >= 80) {
            $use = "YES";
        }
        else {
            $use = "NO";
        }

        #输出结果
        print "$line[1]\t$line[2]\t$tp1\t$tp2\t$tp3\t$bao\t$use\n";
        print OUT "$id\t$y\t$reftaxid\t$amplicon_seq\t$tp1\t$tp2\t$tp3\t$bao\t$use\n";
    }
    close OUT;
    close PRI;
    close NORUN;

    #输出包容性表格
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
}

if ($type eq "Both") {
    if ($runname) {
        open OUT1, ">$outdir/$runname.tNGS.primer.check.xls" or die $!;
        open OUT2, ">$outdir/$runname.qPCR.primer.check.xls" or die $!;
    }
    else {
        open OUT1, ">$outdir/tNGS.primer.check.xls" or die $!;
        open OUT2, ">$outdir/qPCR.primer.check.xls" or die $!;
    }
    #外引物
    open PRI, "<$outdir/External_Primers.csv" or die $!;
    open NORUN, ">$outdir/norun-tNGS-Primers.csv" or die $!; #输出有问题的引物对
    while (<PRI>) {
        chomp;
        #next if $.==1;
        if ($. == 1) {
            my @line = split /,/, $_;
            $line[0] = "ID";
            my $x = join("\t", @line[0 .. 15]);
            print OUT1 "$x\t";
            print OUT1 "Taxid\tAmplicon_Seq\t1-specificity\t2-specificity\tSE751-specificity\tSE752-specificity\tSE1001-specificity\tSE1002-specificity\tinclusiveness\tuse\n";
            print NORUN "$_\n";
        }
        next if $. == 1;
        #print "$_\n";
        my @line = split /,/, $_;
        for (my $a = 0; $a <= $#line; $a++) {
            $line[$a] =~ s/\.0//;
        }
        my $y = join("\t", @line[1 .. 15], @line[18 .. 19]);
        my $primer = join("-", @line[1 .. 2]);
        # my @line2 = split/_/,$line[3];
        # my $species = join(" ",@line2[0..$#line2-2]);
        my @line2 = split /-/, $line[3];
        my $species = $line2[0];
        my $reftaxid;     #物种taxid
        my $id;           #引物唯一ID
        my $amplicon_seq; #扩增子序列
        my ($tp1, $tp2, $tp3, $tp4, $tp5, $tp6, $bao);

        if ($species2taxid{$species}) {
            $reftaxid = $species2taxid{$species};
            print "$species\t$reftaxid\n";
            `echo $primer $reftaxid > id.txt`;
            $id = &md5("id.txt");
            `rm id.txt`;

            #创建物种包含的所有子taxid的对应表
            `taxonkit list --show-rank --show-name --ids $reftaxid --data-dir /home/yehui/software/TaxonKit > idtemp.txt`;
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
                if ($gcf_num <= 2) {
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

            #引物序列1
            open OU1, ">primer1.fa" or die $!;
            print OU1 ">$species\n$line[1]";
            close OU1;
            &blast("primer1.fa", $blastdatabase);
            $tp1 = &parseblast("out.blast", $species, \%taxid2species);

            #引物序列2
            open OU2, ">primer2.fa" or die $!;
            print OU2 ">$species\n$line[2]";
            close OU2;
            &blast("primer2.fa", $blastdatabase);
            $tp2 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序75bp一端
            open OU3, ">SE751.bed" or die $!;
            my $add1 = 75 - $line[12];
            my $SE751_end = $line[5] + $add1;
            print OU3 "$line[4]\t$line[5]\t$SE751_end\n";
            close OU3;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE751.bed");
            &blast("out.fa", $blastdatabase);
            $tp3 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序75bp另一端
            open OU4, ">SE752.bed" or die $!;
            my $add2 = 75 - $line[13];
            my $SE752_start = $line[6] - $add2;
            print OU4 "$line[4]\t$SE752_start\t$line[6]\n";
            close OU4;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE752.bed");
            &blast("out.fa", $blastdatabase);
            $tp4 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序100bp一端
            open OU5, ">SE1001.bed" or die $!;
            my $add3 = 100 - $line[12];
            my $SE1001_end = $line[5] + $add3;
            print OU5 "$line[4]\t$line[5]\t$SE1001_end\n";
            close OU5;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE1001.bed");
            &blast("out.fa", $blastdatabase);
            $tp5 = &parseblast("out.blast", $species, \%taxid2species);

            #SE测序100bp另一端
            open OU6, ">SE1002.bed" or die $!;
            my $add4 = 100 - $line[13];
            my $SE1002_start = $line[6] - $add4;
            print OU6 "$line[4]\t$SE1002_start\t$line[6]\n";
            close OU6;
            &getfa("$allrefdir/$reftaxid/ref.fna", "SE1002.bed");
            &blast("out.fa", $blastdatabase);
            $tp6 = &parseblast("out.blast", $species, \%taxid2species);

            #扩增子序列（包含引物）
            open OU7, ">amp.bed" or die $!;
            print OU7 "$line[4]\t$line[5]\t$line[6]\n";
            close OU7;
            &getfa("$allrefdir/$reftaxid/ref.fna", "amp.bed");
            $amplicon_seq = `sed -n '2p' out.fa`;
            $amplicon_seq =~ s/\n//;

            #留空
            #留空

            #引物包容性分析
            if ($gcf_num >= 10) {
                if (-e "$allrefdir/$reftaxid/refseq.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/refseq.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $species2num{$species} * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "refseq.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
            else {
                if (-e "$allrefdir/$reftaxid/all.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/all.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $gcf_gca_num * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "all.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
        }
        else {
            print "物种拉丁名错误，请检查！\n";
            $reftaxid = "xxx";
            print "$species\t$reftaxid\n";
            print NORUN "$_\n";

            $tp1 = 0;
            $tp2 = 0;
            $tp3 = 0;
            $tp4 = 0;
            $tp5 = 0;
            $tp6 = 0;
            $bao = 0;
        }

        #判断引物对是否可用
        my $use;
        if (($tp1 >= 80 || $tp2 >= 80 || $tp3 >= 80 || $tp4 >= 80 || $tp5 >= 80 || $tp6 >= 80) && $bao >= 80) {
            $use = "YES";
        }
        else {
            $use = "NO";
        }

        #输出结果
        print "$line[1]\t$line[2]\t$tp1\t$tp2\t$tp3\t$tp4\t$tp5\t$tp6\t$bao\t$use\n";
        print OUT1 "$id\t$y\t$reftaxid\t$amplicon_seq\t$tp1\t$tp2\t$tp3\t$tp4\t$tp5\t$tp6\t$bao\t$use\n";
    }
    close OUT1;
    close PRI;
    close NORUN;

    #输出包容性表格
    my $landir1 = "$outdir/2lan/tNGS";
    `mkdir -p $landir1`;
    foreach my $x (keys %taxid2primer) {
        #taxid
        open TAX, ">$landir1/$x.xls" or die $!;
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

    #内引物
    open PRI, "<$outdir/Internal_Primers.csv" or die $!;
    open NORUN, ">$outdir/norun-qPCR-Primers.csv" or die $!; #输出有问题的引物对
    while (<PRI>) {
        chomp;
        #next if $.==1;
        if ($. == 1) {
            my @line = split /,/, $_;
            $line[0] = "ID";
            my $x = join("\t", @line);
            print OUT2 "$x\t";
            print OUT2 "Taxid\tAmplicon_Seq\t1-specificity\t1-specificity\t2-specificity\tamplicon\tinclusiveness\tuse\n";
            print NORUN "$_\n";
        }
        next if $. == 1;
        #print "$_\n";
        my @line = split /,/, $_;
        for (my $a = 0; $a <= $#line; $a++) {
            $line[$a] =~ s/\.0//;
        }
        my $y; #引物设计数据
        my $length = (scalar @line) - 1;
        if ($length == 17) {
            $y = join("\t", @line[1 .. 17]);
        }
        else {
            $y = join("\t", @line[1 .. 16], "-");
        }
        #my $y = join("\t",@line);
        my $primer = join("-", @line[1 .. 2]);
        # my @line2 = split/_/,$line[3];
        # my $species = join(" ",@line2[0..$#line2-2]);
        my @line2 = split /-/, $line[3];
        my $species = $line2[0];
        my $reftaxid;     #物种taxid
        my $id;           #引物唯一ID
        my $amplicon_seq; #扩增子序列
        my ($tp1, $tp2, $tp3, $bao);

        if ($species2taxid{$species}) {
            $reftaxid = $species2taxid{$species};
            print "$species\t$reftaxid\n";
            `echo $primer $reftaxid > id.txt`;
            $id = &md5("id.txt");
            `rm id.txt`;

            #创建物种包含的所有子taxid的对应表
            `taxonkit list --show-rank --show-name --ids $reftaxid --data-dir /home/yehui/software/TaxonKit > idtemp.txt`;
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
                if ($gcf_num <= 2) {
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

            #引物序列1
            open OU1, ">primer1.fa" or die $!;
            print OU1 ">$species\n$line[1]";
            close OU1;
            &blast("primer1.fa", $blastdatabase);
            $tp1 = &parseblast("out.blast", $species, \%taxid2species);

            #引物序列2
            open OU2, ">primer2.fa" or die $!;
            print OU2 ">$species\n$line[2]";
            close OU2;
            &blast("primer2.fa", $blastdatabase);
            $tp2 = &parseblast("out.blast", $species, \%taxid2species);

            #插入片段
            open OU3, ">amplicon.bed" or die $!;
            print OU3 "$line[4]\t$line[8]\t$line[9]\n";
            close OU3;
            &getfa("$allrefdir/$reftaxid/ref.fna", "amplicon.bed");
            &blast("out.fa", $blastdatabase);
            $tp3 = &parseblast("out.blast", $species, \%taxid2species);

            #扩增子序列（包含引物）
            open OU4, ">amp.bed" or die $!;
            print OU4 "$line[4]\t$line[5]\t$line[6]\n";
            close OU4;
            &getfa("$allrefdir/$reftaxid/ref.fna", "amp.bed");
            $amplicon_seq = `sed -n '2p' out.fa`;
            $amplicon_seq =~ s/\n//;

            #引物包容性分析
            if ($gcf_num >= 10) {
                if (-e "$allrefdir/$reftaxid/refseq.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/refseq.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $species2num{$species} * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "refseq.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
            else {
                if (-e "$allrefdir/$reftaxid/all.fna" && %chrom2genom) {
                    &getamplicon("$allrefdir/$reftaxid/all.fna", $line[1], $line[2]);
                    my ($len, $uni) = &parseamplicon("amplicon.txt", \%chrom2genom);
                    $bao = $len / $gcf_gca_num * 100;
                    $bao = sprintf "%.3f", $bao;

                    my @unigenom = @$uni;
                    foreach my $i (@unigenom) {
                        push @{$taxid2primer{$reftaxid}{$primer}}, $i;
                    }
                }
                else {
                    print "all.fna数据库不存在，请检查！\n";
                    $bao = 0;
                }
            }
        }
        else {
            print "物种拉丁名错误，请检查！\n";
            $reftaxid = "xxx";
            print "$species\t$reftaxid\n";
            print NORUN "$_\n";

            $tp1 = 0;
            $tp2 = 0;
            $tp3 = 0;
            $bao = 0;
        }

        #判断引物对是否可用
        my $use;
        if (($tp1 >= 80 || $tp2 >= 80 || $tp3 >= 80) && $bao >= 80) {
            $use = "YES";
        }
        else {
            $use = "NO";
        }

        #输出结果
        print "$line[1]\t$line[2]\t$tp1\t$tp2\t$tp3\t$bao\t$use\n";
        print OUT2 "$id\t$y\t$reftaxid\t$amplicon_seq\t$tp1\t$tp2\t$tp3\t$bao\t$use\n";
    }
    close OUT2;
    close PRI;
    close NORUN;

}

#=======================子函数==========================#
sub blast {
    #$_[0]：需要blast的fa序列
    #$_[1]：blast用到的数据库
    #`blastn -task blastn-short -query $_[0] -db $_[1] -num_threads 50 -word_size 7 -outfmt "7 staxid ssciname qseqid sseqid length qlen slen sstart send qstart qend pident nident evalue bitscore" -out out.blast`;
    `blastn -task blastn-short -query $_[0] -db $_[1] -num_threads 100 -word_size 7 -outfmt "6 staxid ssciname qlen mismatch qstart qend pident nident evalue bitscore" -out out.blast`;
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
        print "blast无结果\t";
        print "$hits\n";
        $p = "0";
    }

    $p = sprintf "%.3f", $p;
    print "$p\n";
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
    `/home/yehui/software/bin/seqkit amplicon $_[0] -F $_[1]  -R $_[2]  -m 2 -j 20 --bed >amplicon.txt`; #设置错误匹配2
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
