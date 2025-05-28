#!/usr/bin/perl -w 
#########################################################
# Author: stone
# Created Time : 2023.05.25 10:00:00 AM CST
# File Name: primer_design.pl
# Version: v1.0
#########################################################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::MoreUtils ':all';
use Cwd qw(getcwd abs_path);


# 帮助文档
=head1 Description

    合并内引物和外引物，将外引物进行标记，将外引物插入片段反向，并输出结果表格。

=head1 Usage

    $0 -n <物种名> -i <内引物文件> -o <外引物文件> -r <参考基因组>

=head1 Parameters
    -n  [str]   物种中文名 
    -i  [str]   内引物文件
    -o  [str]   外引物文件
    -r  [str]   参考基因组
    -l  [int]   左右延长碱基长度
    -h  [str]   帮助信息
=cut

my ($name,$inside,$outside,$ref,$len,$outdir,$help);
GetOptions(
    "n|name=s"      => \$name,
    "i|inside=s"    => \$inside,
    "o|outside=s"   => \$outside,
    "r|refseq=s"    => \$ref,
    "l|length:i"    => \$len,
    "d|outdir:s"    => \$outdir,
    "h|help:s"      => \$help
    );

die `pod2text $0` if ((!$inside && !$outside) or ($help));

$len ||= 20;

$outdir ||= "./";
`mkdir -p $outdir`;
$outdir = abs_path($outdir);

# 读取内引物信息
open FH1, "<$inside" or die $!;
my $header1;
my %table1;
while(<FH1>){
    chomp;
    my @line = split/\t/,$_;
    if($.==1){
        $header1 = $_;
    }
    next if $.==1;
    if($_ =~ /YES/){
        $table1{$line[4]} = $_;
    }else{
        next;
    }
}
close FH1;

# 读取外引物信息
open FH2, "<$outside" or die $!;
my $header2;
my %table2;
while(<FH2>){
    chomp;
    my @line = split/\t/,$_;
    if($.==1){
        $header2 = $_;
    }
    next if $.==1;
    if($_ =~ /YES/){
        my $len_F = length($line[1]);#F引物长度
        my $len_R = length($line[2]);#R引物长度
        my $len_B = $line[7] - $len_F - $len_R;#除去引物后长度

        my $B_start = $len_F;#中间片段开始位置
        my $B_end = $len_F + $len_B;#中间片段结束位置

        my $F = substr($line[19],0,$len_F);
        my $Block = substr($line[19],$B_start,$len_B); #中间片段序列
        my $R = substr($line[19],$B_end,$len_R);

        my $Block_rev = reverse($Block);#反向中间片段序列
        my $amp = join("","[",$line[1],"]",$Block,"[",$R,"]");#扩增子标记后序列
        $line[19] = $amp;
        $table2{$line[3]} = join("\t",@line,$Block_rev);
    }else{
        next;
    }
}
close FH2;

# print Dumper (%table2);
# print Dumper (%table1);

# 输出连接后的表格
# 遍历 table2，将其中的行与 table1 中的行连接
open OUT,">$outdir/$name\_merge_primer_tmp.xls" or die $!;
print OUT "species\t$header1\t$header2\tBlock_rev\n";
foreach my $i(sort keys %table2){
    if(exists $table1{$i}){
        print OUT "$name\t$table1{$i}\t$table2{$i}\n";
    }else{
        next;
    }
}
close OUT;

open IN,"<$outdir/$name\_merge_primer_tmp.xls" or die $!;
open OUT,">$outdir/$name\_merge_primer.xls" or die $!;
while(<IN>){
    chomp;
    my @line = split/\t/,$_;
    if($.==1){
        my $x = join("\t",@line[0..19],"Amplicon_Seq_add","Amplicon_Seq","in_primer_pair_specificity","in_primer_specificity","in_primer_inclusiveness","probe-specificity","probe_inclusiveness",@line[34..53],"Amplicon_Seq_add_before","Amplicon_Seq","Amplicon_Seq_add_after","out_primer_pair_specificity","out_primer_specificity","out_primer_inclusiveness",@line[64,65],"level");
        print OUT "$x\n";
    }else{
        #内引物1base+1base
        #内引物扩增片段左端延长
        my $a = $line[7] - 1 - $len;
        #内引物扩增片段右端延长
        my $b = $line[8] + $len;
        #生成bed文件,获取序列
        open QPCR,">qpcr_add.bed" or die $!;
        print QPCR "$line[6]\t$a\t$b\n";
        close QPCR;
        &getfa($ref,"qpcr_add.bed");
        my $Amplicon_Seq_add = `sed -n '2p' out.fa`;
        $Amplicon_Seq_add =~ s/\n//;
        `rm qpcr_add.bed out.fa`;

        #外引物0base+0base
        ### 
        #外引物扩增片段左端延长
        # my $c = $line[40] - $len ;
        # my $d = $line[40];
        my $c = $line[40] - $len - 1;
        my $d = $line[40] - 1;
        #生成bed文件,获取序列
        open BF,">tngs_add_before.bed" or die $!;
        print BF "$line[39]\t$c\t$d\n";
        close BF;
        &getfa($ref,"tngs_add_before.bed");
        my $Amplicon_Seq_add_before = `sed -n '2p' out.fa`;
        $Amplicon_Seq_add_before =~ s/\n//;
        #外引物扩增片段右端延长
        # my $e = $line[41] + 1;
        # my $f = $line[41] + $len + 1;
        my $e = $line[41];
        my $f = $line[41] + $len;
        #生成bed文件,获取序列
        open AF,">tngs_add_after.bed" or die $!;
        print AF "$line[39]\t$e\t$f\n";
        close AF;
        &getfa($ref,"tngs_add_after.bed");
        my $Amplicon_Seq_add_after = `sed -n '2p' out.fa`;
        $Amplicon_Seq_add_after =~ s/\n//;
        `rm tngs_add_before.bed tngs_add_after.bed out.fa`;

        my ($inprimer_pair_primer_spe,$in_f_per,$in_r_per,$inprimer_probe_spe,$inprimer_primer_inc,$inprimer_probe_inc) = @line[21,23,25,27,29,31];
        my $inprimer_primer_spe = "$in_f_per|$in_r_per";
        my $in_evaluate = join("\t",$inprimer_pair_primer_spe,$inprimer_primer_spe,$inprimer_primer_inc,$inprimer_probe_spe,$inprimer_probe_inc);

        my ($outprimer_primer_pair_spe,$out_f_per,$out_r_per,$outprimer_primer_inc) = @line[55,57,59,61];
        my $outprimer_primer_spe = "$out_f_per|$out_r_per";
        my $out_evaluate = join("\t",$outprimer_primer_pair_spe,$outprimer_primer_spe,$outprimer_primer_inc);

        my $level = "A";
        if($inprimer_pair_primer_spe >= 95 && $inprimer_primer_inc >= 95 && $inprimer_probe_inc >= 95){
            $level = "S";
        }elsif($inprimer_pair_primer_spe >= 90 && $inprimer_primer_inc >= 90 && $inprimer_probe_inc >= 90){
            $level = "A";
        }elsif($inprimer_pair_primer_spe >= 80 && $inprimer_primer_inc >= 80 && $inprimer_probe_inc >= 80){
            $level = "B";
        }elsif($inprimer_pair_primer_spe >= 60 && $inprimer_primer_inc >= 60 && $inprimer_probe_inc >= 60){
            $level = "C";
        }else{
            $level = "O";
        }
        my $inprimer_seq = $line[20];
        my $outprimer_seq = $line[54];

        my $y = join("\t",@line[0..19],$Amplicon_Seq_add,$inprimer_seq,$in_evaluate,@line[34..53],$Amplicon_Seq_add_before,$outprimer_seq,$Amplicon_Seq_add_after,$out_evaluate,@line[64,65],$level);
        print OUT "$y\n";
    }
}
# `rm $name\_merge_primer_tmp.xls`;

sub getfa{
    #$_[0]：参考基因组
    #$_[1]：bed区间文件
    `/home/chenwenjing/pipeline/IDseq_v2.0/bin/bedtools getfasta -fi $_[0] -bed $_[1] -fo out.fa`;
}