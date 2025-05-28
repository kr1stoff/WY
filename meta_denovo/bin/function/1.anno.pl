#!/usr/bin/perl -w
use strict;

die "perl $0 <blast_m8_file> <output_file> <db_info_file>" if (@ARGV != 3);

my ($infile, $outfile, $info_file) = @ARGV;

open IN,"$info_file" or die;
my %info = ();
while(my $line = <IN>){
	chomp $line;
	my @array = split(/\s+/,$line, 2);
	$info{$array[0]} = $array[1];
};
close IN;

open IN, "$infile" or die;
open OUT, ">$outfile" or die;
print OUT "Gene_id\tIdentity\tE_value\tSubject_id\tSubject_description\n";
while(my $line = <IN>){
	chomp $line;
	next if ($line =~ /^Gene_id/);
	my ($que,$sub,$id,$eval) = (split(/\t/,$line))[0,1,2,10];
	my $tmp = "";
	if(!exists $info{$sub}){
		$tmp = "NA";
	}else{
		$tmp = $info{$sub};
	};
	print OUT "$que\t$id\t$eval\t$sub\t$tmp\n";
};
close IN;
close OUT;

