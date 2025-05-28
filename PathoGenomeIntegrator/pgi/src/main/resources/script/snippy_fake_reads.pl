#!/usr/bin/env perl

# @Author   : mengxf
# @Version  : 1.0
# @Create   : 2024/03/22
# @Modified : 2024/04/03
# https://github.com/tseemann/snippy/blob/master/bin/snippy

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename;

# 变量定义
my $in_fasta;
my $out_fastq1;
my $out_fastq2;
my $min_cov    = 100;    # 默认最小覆盖
my $read_len   = 500;    # 默认读取长度
my $paired_end = 0;      # 默认单端

# 参数解析和验证
GetOptions(
    'infa=s'    => \$in_fasta,
    'outfq1=s'  => \$out_fastq1,
    'outfq2=s'  => \$out_fastq2,
    'mincov=i'  => \$min_cov,
    'readlen=i' => \$read_len,
    'help'      => sub { HelpMessage(0) },
) or HelpMessage(1);

if ( !defined($in_fasta) || !defined($out_fastq1) ) {
    HelpMessage( 1,
        "must specify the input FASTA file and the output read1 FASTQ file." );
}

# 时间戳,文件名 变量
my $timestamp = time;
my $basen     = basename($in_fasta);
$basen =~ s/\./_/g;

# 主流程
if ( !defined($out_fastq2) ) {
    my $paired_end = 0;
}
else {
    my $paired_end = 1;
}

generate_fake_reads(
    $in_fasta, $out_fastq1, $out_fastq2, $min_cov,
    $read_len, $basen,      $timestamp,  $paired_end
);

# 错误消息改进
sub HelpMessage {
    my ( $exit_code, $message ) = @_;
    if ($message) {
        print "ERROR: $message\n\n";
    }
    print <<"END_HELP";
Usage: $0 [options]

Options:
    --infa=FILE         the path to the input FASTA file.
    --outfq1=FILE       the path to the output READ1 FASTQ file.
    --outfq2=FILE       the path to the output READ2 FASTQ file, if provided, output in PE mode.
    --mincov=INT        the minimum coverage of the fake reads. (default: 100)
    --readlen=INT       the length of the fake reads. (default: 500)
    --help              show this help message.
END_HELP
    exit($exit_code);
}

sub generate_fake_reads {
    my (
        $in_fasta, $out_fastq1, $out_fastq2, $min_cov,
        $read_len, $basen,      $timestamp,  $paired_end
    ) = @_;
    my $fake_read_cov = 2 * $min_cov;
    my $stride        = $read_len / ( 0.5 * $fake_read_cov );
    my @out_handles;

    if ($paired_end) {
        open $out_handles[0], '>:encoding(UTF-8)', $out_fastq1
          or die "Cannot open $out_fastq1: $!";
        open $out_handles[1], '>:encoding(UTF-8)', $out_fastq2
          or die "Cannot open $out_fastq2: $!";
    }
    else {
        open $out_handles[0], '>', $out_fastq1
          or die "Cannot open $out_fastq1: $!";
    }

    my $in = Bio::SeqIO->new( -file => $in_fasta, -format => 'fasta' );

    while ( my $seq = $in->next_seq ) {
        my @dna    = ( uc( $seq->seq ), uc( $seq->revcom->seq ) );
        my $L      = $seq->length;
        my $len    = ( $L < $read_len ) ? $L : $read_len;
        my $counter = 0;

        for (
            my $start_pos = -$len ;
            $start_pos < $L + $len ;
            $start_pos += $stride
          )
        {
            my $start = ( $start_pos > 0 ) ? $start_pos : 0;
            $start = ( $start < ( $L - $len ) ) ? $start : ( $L - $len );
            $counter++;

            if ($paired_end) {
                print { $out_handles[0] }
                  "\@$basen:$timestamp:read${counter} 1/2\n",
                  substr( $dna[0], int($start), $len ), "\n", "+\n",
                  ('F') x $len,
                  "\n";
                print { $out_handles[1] }
                  "\@$basen:$timestamp:read${counter} 2/2\n",
                  substr( $dna[1], int($start), $len ), "\n", "+\n",
                  ('F') x $len,
                  "\n";
            }
            else {
                print { $out_handles[0] } "\@$basen:$timestamp:read$counter\n",
                  substr( $dna[0], int($start), $len ), "\n", "+\n",
                  ('F') x $len,
                  "\n";
            }
        }
    }
    close $_ for @out_handles;    # Close all file handles after use
}
