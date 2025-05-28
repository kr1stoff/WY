use strict;

my $blast = shift;

open A,$blast || die $!;
while(<A>){
	chomp;
	my @a = split /\t/,$_;
	if ($a[2]>=90 && $a[14]>=80){
		print "$_\n";
	}
}
close A;

