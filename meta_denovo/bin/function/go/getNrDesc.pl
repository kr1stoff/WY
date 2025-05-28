#!/sur/bin/perl -w

=pod
description: get desc from nr
author: Zhang Fangxian, zhangfx@genomics.cn
created date: 20090915
modified date: 20101127, 20100612, 20091105, 20091029, 20090916
=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);

my ($input, $nr, $rank, $output, $help);

GetOptions("input:s" => \$input, "nr:s" => \$nr, "rank:i" => \$rank, "output:s" => \$output, "help|?" => \$help);
$rank = 1 if (!defined $rank);
if (!defined $input || !defined $nr || defined $help) {
	die << "USAGE";
description: get desc from nr
usage: perl $0 [options]
options:
	-input <file> *  input file, blastout of against nr in m8 format
	-nr <file> *     nr database
	-rank <int>      rank cutoff for valid hit from blastall, default is 1
	-output <file>   output tabular file, default is ./[-input].xls

	-help|?          help information
e.g.:
	perl $0 -input 1.tab -nr $nr -rank 5 -output 1.xls
USAGE
}

if (!-f $input) {
	die "file $input not exists\n";
}

$output ||= &getFileName($input) . ".xls";

my (%result, %nrs);
open IN, "< $input" or die $!;
while (<IN>) {
	my @tabs = split /\t/, $_;
	$tabs[0] = (split /\s/, $tabs[0])[0];
	if (exists $result{$tabs[0]} && @{$result{$tabs[0]}} < $rank) {
		push @{$result{$tabs[0]}}, [$tabs[1], $tabs[10]];
		$nrs{$tabs[1]} = 1;
	} elsif (!exists $result{$tabs[0]}) {
		push @{$result{$tabs[0]}}, [$tabs[1], $tabs[10]];
		$nrs{$tabs[1]} = 1;
	}
}
close IN;

open NR, "< $nr" or die $!;
while (<NR>) {
	s/[\r\n]//g;
	next unless (/^>/);
	$_ =~ s/^>//;
	for my $i (split /\s*>/, $_) {
	my @tabs = split /\s/, $i, 2;
		if (exists $nrs{$tabs[0]}) {
			$nrs{$tabs[0]} = $tabs[1];
		}
	}
}
close NR;

open OUT, "> $output" or die $!;
for my $gene (keys %result) {
	for (0 .. $#{$result{$gene}}) {
		last if ($nrs{$result{$gene}->[0]->[0]} !~ /predicted|hypothetical/i);
		my $temp = shift @{$result{$gene}};
		push @{$result{$gene}}, $temp;
	}
	print OUT "$gene\t$result{$gene}->[0]->[0]\t$result{$gene}->[0]->[1]\t$nrs{$result{$gene}->[0]->[0]}\n";
}
close OUT;

exit 0;

sub getFileName {
	my ($file_name) = @_;
	$file_name = (split /[\/\\]/, $file_name)[-1];
	$file_name =~ s/\.[^\.]*$//;
	return $file_name;
}
