#!/usr/bin/perl
use warnings;
use threads;
use FindBin qw($Bin $RealBin);
use JSON;
use File::Basename qw(basename);
use YAML qw(LoadFile);
use Getopt::Long;

# 读取所用软件
my $f_software = "$RealBin/../config/software.yml";
our $app = LoadFile($f_software);
my $fastp = $app->{"fastp"};

# 命令行参数
my %opts;
GetOptions(\%opts, "r1=s", "r2=s", "o=s", "l=s", "q=s", "p=s", "min=i", "tf=i", "tr=i");

$usage = <<"USAGE";
        Usage:  perl $0 -r1 <read_1.fq.gz> -r2 <read_2.fq.gz> [options]
        Notes:  !!! quality 33,q20,50%; N 10%; shortest 50bp; remove adapater; threads 4 !!!
        Options:
                -r1     strings     read1 path
                -r2     strings     read2 path
                -l      strings     read length, default 150
                -o      prefix      output prefix
                -q      quality     [33]/64
                -p      program     fastp path
                -min    int         min length, reads shorter than length_required will be discarded
                -tf     int         length of forward primer trim
                -tr     int         length of reverse primer trim
        eg:
                perl $0 -r1 read_1.fq.gz -r2 read_2.fq.gz -o /xx/oo/prefix <...>
USAGE

# Main
die $usage if (!$opts{r1} or !$opts{o});
$opts{l} = $opts{l} ? $opts{l} : 150;
$opts{q} //= 33;
$opts{min} //= 50;
my $N = int($opts{l} * 0.1);
my $phred = $opts{q} == 64 ? '--phred64' : '';
$opts{other_opt} = scalar(@ARGV) > 0 ? join(" ", @ARGV) : '';
if (exists $opts{tf} || exists $opts{tr}) {
    $opts{tf} //= 0;
    $opts{tr} //= 0;
    $opts{other_opt} = "-f $opts{tf} -F $opts{tr} $opts{other_opt}";
}

$fastp = exists $opts{p} ? $opts{p} : $fastp;
my $adapter = 'AGATCGGAAGAGC';

my $cmd_part = "$fastp -a $adapter -q 20 -u 50 -n $N -l $opts{min} -w 1 $phred -j $opts{o}.json -h $opts{o}.html $opts{other_opt}";
if (exists $opts{r2}) {
    Excu("$cmd_part -i $opts{r1} -I $opts{r2} -o $opts{o}_1.fq.gz -O $opts{o}_2.fq.gz");
}
else {
    Excu("$cmd_part -i $opts{r1} -o $opts{o}_1.fq.gz");
}

## get stat from json
open IN, "<$opts{o}.json" or die "Can't open the file! $!";
while (<IN>) {$json_text .= $_;}
close IN;
my $json = decode_json $json_text;
my $before = $json->{summary}->{before_filtering};
my $after = $json->{summary}->{after_filtering};
my $sample = basename($opts{o});
my $be_gc = $before->{gc_content} * 100;
my $be_q20 = $before->{q20_rate} * 100;
my $be_q30 = $before->{q30_rate} * 100;
my $af_gc = $after->{gc_content} * 100;
my $af_q20 = $after->{q20_rate} * 100;
my $af_q30 = $after->{q30_rate} * 100;
my $lowq = $before->{total_reads} - $after->{total_reads};

open OUT, ">$opts{o}.stat.xls" or die "Can't write the file! $!";
print OUT <<STAT;
Sample	fastq_file	read_length	raw_GC_rate	raw_Q20_rate	raw_reads	raw_bases	adapter_reads	low_qual_reads	clean_GC_rate	clean_Q20_rate	clean_reads	clean_bases
$sample	$sample	$before->{read1_mean_length}	$be_gc	$be_q20	$before->{total_reads}	$before->{total_bases}	0	$lowq	$af_gc	$af_q20	$after->{total_reads}	$after->{total_bases}
STAT
close OUT;

#### Some Function
sub Excu {
    my $cmd = shift;
    print "$cmd\n";
    system($cmd);
    return;
}
