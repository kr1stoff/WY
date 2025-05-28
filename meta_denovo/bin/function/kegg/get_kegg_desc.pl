use strict;
use Data::Dumper;

my $desc = shift;
my $eggNOG2KEGG = shift;

my %hash;
open A,$desc || die $!;
while(<A>){
	chomp;
	my @a = split /\t/,$_;
	$hash{$a[0]} = $a[1];
}
close A;

# GMKBENJP_00001	ko:K02470,ko:K02622

print "基因\tKO\t描述\n";

open B,$eggNOG2KEGG || die $!;
while(<B>){
	chomp;
	next if /^#/;
	my @b = split;
	my @c = split/,/,$b[1];
	my @desc;
	foreach my $ko (@c){
		my $k = $1 if $ko =~ /ko:(.*)/;
		my $desc = $hash{$k};
		push @desc,$desc;
	}
	print "$b[0]\t$b[1]\t";
	print join (", ",@desc)."\n";
}
