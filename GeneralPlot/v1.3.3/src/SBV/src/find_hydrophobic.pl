#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: find_hydrophobic.pl
#
#        USAGE: ./find_hydrophobic.pl  
#
#  DESCRIPTION: search the hydrophobic aa in protein sequence and extract the 
#               position. 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Peng Ai (), apanhui@gmail.com
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 02/21/2014 11:04:38 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

use Bio::SeqIO;

my $hydrophobic_aa = "YWFVLIPA";
my @array = split // , $hydrophobic_aa;
my %hash = map {$_ => 1} @array; 

my $file = shift @ARGV;
my $region = shift @ARGV;

my $inseq = Bio::SeqIO->new(-file=>$file,-format=>"Fasta");
my $seq = $inseq->next_seq;

my $name = $seq->id;
my $seqstr = $seq->seq;
my @aas = split // , $seqstr;
my @res = map { $hash{$_} ? 1 : 0 } @aas;
my $str = join "" , @res;

my @region = read_region($region);
for (@region)
{
	my ($sta,$end) = @$_;
	my $len = $end - $sta + 1;
	my $target = substr($str,$sta-1,$len);
	while($target =~ /1+/g)
	{
		my $hy_sta = pos($target) - length($&) + $sta;
		my $hy_end = pos($target) + $sta - 1;
		print "$hy_sta\t$hy_end\n";
	}
}

sub read_region
{
	my $file = shift;
	my @region;

	open FH,$file or die;
	while(<FH>)
	{
		chomp;
		my ($sta,$end) = split;
		($sta,$end) = ($end,$sta) if ($sta > $end);
		push @region , [$sta,$end];
	}
	close FH;

	return @region;
}
