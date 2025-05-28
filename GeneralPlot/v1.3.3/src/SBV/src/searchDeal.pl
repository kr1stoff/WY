#!/usr/bin/perl -w
# Name: searchDeal.pl
# 
# Func: used to deal search result file, like BLAST/FASTA etc.
# 1. format transform
# 2. result filt
# 3. result visualization
#
# Author: aipeng
# E-mail: ap.anhui@gmail.com
#
# Version
# ===========================
# 1.0 beta  DATE:12-26-2012


use strict;

use Bio::SearchIO;
use Getopt::Std;
use Math::Round;

my %opts = (f=>'blast',o=>'sbv',k=>'KEYNAME',D=>'.',E=>10,S=>0,I=>0,L=>0,C=>0,c=>0);
getopts('f:o:BC:c:L:I:E:S:k:D:',\%opts);

&usage if (@ARGV == 0);

my $infile = shift @ARGV;
my %length;
if ($opts{f} eq "blasttable"){readLst(shift @ARGV,\%length)}

## creat output directory and output file handle
mkdir $opts{D} unless -d $opts{D};
my $ofile = "$opts{D}/$opts{k}.$opts{o}";
my $ofh;
open $ofh , ">" , $ofile or die "can't open file: $ofile";

## read input file
# init the searchio par 
my %par = (-format=>$opts{f},-file=>$infile,-check_all_hits=>1);
$par{'-best_hit_only'} = 1 if ($opts{B});
my $in = new Bio::SearchIO(%par);

while(my$result = $in->next_result)
{
	my $query_id  = $result->query_name();
    my $query_len = $opts{f} eq "blasttable" ? $length{$query_id} : $result->query_length();
	my $query_desc = $result->query_description || "--";

    while(my$hit = $result->next_hit)
    {
		my $sbj_id = $hit->name;
		my $sbj_len = $opts{f} eq "blasttable" ? $length{$sbj_id} : $hit->length;
		my $sbj_desc = $hit->description || "--";
		
		while(my$hsp = $hit->next_hsp)
        {
			my $evalue = $hsp->evalue;
			my $bits = $hsp->bits;
			my $identity = nearest 0.01 , $hsp->percent_identity;
			my $matchLen = $hsp->length('total');
			my $qry_start = $hsp->start('query');
			my $qry_end = $hsp->end('query');
			my $sbj_start = $hsp->start('hit');
			my $sbj_end = $hsp->end('hit');
			my $sbj_strand = $hsp->strand('subject');

			my $scov = nearest 0.01 , 100 * ($sbj_end-$sbj_start+1) / $sbj_len;
			my $qcov = nearest 0.01 , 100 * ($qry_end-$qry_start+1) / $query_len;
			if ($evalue < $opts{E} && $bits > $opts{S} && $identity > $opts{I} && 
				$matchLen > $opts{L} && $scov >= $opts{C} && $qcov >= $opts{c})
			{
				($sbj_start,$sbj_end) = ($sbj_end,$sbj_start) if (-1 == $sbj_strand);
				my @array = ($query_id,$query_len,$qry_start,$qry_end,$sbj_id,$sbj_len,$sbj_start,$sbj_end,
					$matchLen,$scov,$identity,$evalue,$bits,$query_desc,$sbj_desc);
				save($ofh,$opts{o},\@array);
			}
			
			last if ($opts{B});
        }
		
		last if ($opts{B});
    }
}

close $ofh;
#=====================================
# sub program
#=====================================
# read length list file 
sub readLst
{
	my ($file,$hash) = @_;
	open FH,$file or die "can't open file: $file $!";
	while(<FH>)
	{
		chomp;
		my ($name,$len) = split;
		$hash->{$name} = $len;
	}
	close FH;

	return $hash;
}

sub save
{
	my ($ofh,$format,$array) = @_;

	my $record;
	if ($format eq "sbv")
	{
		my @temp = @$array[0,1,2,3,4,5,6,7,10];
		$record = join "\t" , @temp;
	}
	elsif ($format eq "full")
	{
		$record = join "\t" , @$array;
	}
	
	print $ofh "$record\n";
}

# usage direction for help
sub usage
{
    print qq(
Usage:   
         perl searchDeal.pl [options] <search result file> [length list file]

Options: 
         -f STR    input file format, [blast]
         -o STR    output file format, [sbv]

         -k STR    output file keyname, [KEYNAME]
         -D STR    output Directory , [.]

         -B        fetch the besthit only [not]
         -C FLOAT  the minmum match covergae of subject 0~100, [0]
         -c FLOAT  the minmum match covergae of query 0~100, [0]
         -L INT    the minmum match Length vaule 0~length of the seq, [0]
         -I INT    the minmum percentage identity value 0~100, [0]
         -E FLOAT  the maximum evalue, >=0, [10]
         -S INT    the minmum socre value, [0]

Note: 
         the format of input file is defined in the Bio::SearchIO moudle

\n);
    exit 1;
}
