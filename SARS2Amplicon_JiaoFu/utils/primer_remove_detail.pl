#!/usr/bin/perl
use File::Basename;
use FindBin qw($Bin);
use Getopt::Long;

BEGIN{
    push(@INC,$Bin);
};
my $bin = dirname $Bin;

my $usage="usage:
    perl $0 [Options]
Options:
  - Sample name:
    -i   [s] input bam file, sort by read name (samtools sort -n) <required>
    -pr  [s] primer file, five rows, 'chr left-st left-ed right-st right-ed' <required>
    -ml  [s] min length for output <default: 20>
    -mq  [s] min map quality <default: 5>
    -thes   [s] length extend for cut <default: 3>
    -o  [s] outdir <default: ./>
    -s  [s] sample name <default: sample1>
    -samtools [s] samtools bin path <default: $bin/software/samtools>
author:
    Furion
";

unless(@ARGV>0){
    print "$usage";
    exit;
}

GetOptions(
    "i=s" => \$bamfile, # name sort bam
    "pr=s" => \$primerfile,
    "ml=s" => \$minlength,
    "mq=s" => \$minmapquality,
    "thes=s" => \$lengththesold,
    "o=s" => \$outdir,
    "s=s" => \$sample,
    "samtools=s" => \$samtools,
);

# default
unless(defined $minlength){$minlength=20;}
unless(defined $minmapquality){$minmapquality=5;}
$samtools||="$bin/software/samtools";
$lengththesold||=3;
$outdir||="./";
$sample||="sample1";
# read primer 
my %hashp1;
my %hashp2;
if(-e "$primerfile"){
    open IN,"$primerfile";
    while(<IN>){
        chomp;
        unless($_){next;}
        if(m/^\#/){next;}
        my @all = split(/\s+/,$_);
        $all[1]=$all[1]+1-$lengththesold;
        $all[3]=$all[3]+1;
        $all[4]=$all[4]+$lengththesold;
        my $pos = "$all[1]\t$all[2]";
        $hashp1{$all[0]}{$pos}=1;
        $pos = "$all[3]\t$all[4]";
        $hashp2{$all[0]}{$pos}=1;
    }
    close IN;
}
# step1 #
my %hashskip;
my ($reads,$bases,$minbase,$maxbase) = (0,0,500,0);
open IN,"$samtools view -h -F 256 $bamfile |";
open OUT1,">$outdir/$sample.rmprimer_R1.fastq" or die $!;
open OUT2,">$outdir/$sample.rmprimer_R2.fastq" or die $!;
while(my $r1=<IN>){
    chomp $r1;
    if($r1=~/^\@/){next;}
    my $r2=<IN>;chomp $r2;
    my @a = split(/\t/,$r1);
    my @b = split(/\t/,$r2);
    unless($a[0] eq $b[0]){die "Error: $r1 and $r2 id not same\n";}
    # cigar R1 R2 ,plus soft clip
    if($a[8]>0 && $b[8]<0){
        # ana r1
        my $softa = 0;
        if($a[5]=~/^(\d+)S/){$softa=$1;}
        my $lena = length $a[9];
        foreach my $key(keys %{$hashp1{$a[2]}}){
            my ($start,$end) = (split /\t/,$key)[0,1];
            if($a[3]>=$start && $a[3]<=$end){
                $reads++;
                my $l = $end - $a[3] + 1 + $softa;
                if($lena<$l+$minlength){$hashskip{$a[0]} = 1;}
                else{
                    $a[9] = substr($a[9],$l);
                    $a[10] = substr($a[10],$l);
                    $bases = $bases+$l;
                    minmax($l);
                }
                last;
            }
        }
        # ana r2
        my $softb=0;
        if($b[5]=~/(\d+)S$/){$softb=$1;}
        my $lenb = length $b[9];
        $lenb = $lenb - $softb;
        # cal readed
        my ($ins,$del) = (0,0);
        while($b[5]=~/(\d+)I/mg){$ins = $ins + $1;}
        while($b[5]=~/(\d+)D/mg){$del = $del + $1;}
        my $ed = $b[3]+$lenb-$ins+$del-1;
        foreach my $key(keys %{$hashp2{$b[2]}}){
            my ($start,$end) = (split /\t/,$key)[0,1];
            if($ed>=$start && $ed<=$end){
                $reads++;
                my $rl = $ed - $start + 1;
                my $l = $lenb - $rl ;
                if($l<$minlength){$hashskip{$b[0]} = 1;}
                else{
                    $b[9] = substr($b[9],0,$l);
                    $b[10] = substr($b[10],0,$l);
                    $bases = $bases+$rl;
                    minmax($rl);
                }
                last;
            }
        }
        # rever
        $b[9] = rever_hubu($b[9]);
        $b[10] = rever($b[10]);
    }elsif($a[8<0] && $b[8]>0){
        # ana r1
        my $softa=0;
        if($a[5]=~/(\d+)S$/){$softa=$1;}
        my $lena = length $a[9];
        $lena = $lena - $softa;
        # cal readed
        my ($ins,$del) = (0,0);
        while($a[5]=~/(\d+)I/mg){$ins = $ins + $1;}
        while($a[5]=~/(\d+)D/mg){$del = $del + $1;}
        my $ed = $a[3]+$lena-$ins+$del-1;
        foreach my $key(keys %{$hashp2{$a[2]}}){
            my ($start,$end) = (split /\t/,$key)[0,1];
            if($ed>=$start && $ed<=$end){
                $reads++;
                my $rl = $ed - $start + 1;
                my $l = $lena - $rl ;
                if($l<$minlength){$hashskip{$a[0]} = 1;}
                else{
                    $a[9] = substr($a[9],0,$l);
                    $a[10] = substr($a[10],0,$l);
                    $bases = $bases+$rl;
                    minmax($rl);
                }
                last;
            }
        }
        # ana r2
        my $softb = 0;
        if($b[5]=~/^(\d+)S/){$softb=$1;}
        my $lenb = length $b[9];
        foreach my $key(keys %{$hashp1{$b[2]}}){
            my ($start,$end) = (split /\t/,$key)[0,1];
            if($b[3]>=$start && $b[3]<=$end){
                $reads++;
                my $l = $end - $b[3] + 1 + $softb;
                if($lenb<$l+$minlength){$hashskip{$b[0]} = 1;}
                else{
                    $b[9] = substr($b[9],$l);
                    $b[10] = substr($b[10],$l);
                    $bases = $bases+$l;
                    minmax($l);
                }
                last;
            }
        }
        # rever
        $a[9] = rever_hubu($a[9]);
        $a[10] = rever($a[10]);
    }
    if(exists $hashskip{$a[0]} || exists $hashskip{$b[0]}){next;}
    print OUT1 "\@$a[0]\n$a[9]\n+\n$a[10]\n";
    print OUT2 "\@$b[0]\n$b[9]\n+\n$b[10]\n";
}
close IN;
close OUT1;
close OUT2;

print STDERR "Modify reads: $reads\n";
print STDERR "Modify bases: $bases\n";
print STDERR "Min base: $minbase\n";
print STDERR "Max base: $maxbase\n";
# hanshu
sub minmax
{
    my $v = shift;
    if($v<$minbase){$minbase=$v;}
    if($v>$maxbase){$maxbase=$v;}
}

sub rever_hubu
{
    my $seq_s = shift;
    $seq_s =~ tr/atucgACGUT/TAAGCTGCAA/;
    $seq_a = reverse($seq_s);
    return $seq_a;
}

sub rever
{
    my $seq_s = shift;
    $seq_a = reverse($seq_s);
    return $seq_a;
}