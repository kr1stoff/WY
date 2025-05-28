use strict;

my $in = shift;
open A,$in || die $!;
while(<A>){
	chomp;
	next if /^#/;
	my @a = split;
	if ($a[1] eq ""){
		next;
	}
	else{
		my @b = split /,/,$a[1];
		foreach my $go (@b){
			print "$a[0]\t$go\n";
		}
	}

}
