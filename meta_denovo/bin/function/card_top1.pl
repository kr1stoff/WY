use strict;

my $in = shift;
my $num = 1;

my %hash;
open A,$in || die $!;
while(<A>){
	chomp;
	next if(/^比对基因ID/);
	my @a = split;
	$hash{$a[0]} .= $_."__";
}
close A;

foreach my $query (sort keys %hash){
	my @tmp = split /__/,$hash{$query};
	if ($#tmp<$num-1){
		for (my $i=0;$i<=$#tmp;$i++){
			print "$tmp[$i]\n";
		}
	}
	else{
		for (my $i=0;$i<=$num-1;$i++){
			print "$tmp[$i]\n";
		}
	}
}
