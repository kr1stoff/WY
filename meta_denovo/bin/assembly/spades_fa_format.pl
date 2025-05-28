use strict;

my $fa = shift;

## >NODE_1_length_337377_cov_239.997910

#my $id_old;
open A,$fa ||die $!;
while(<A>){
	chomp;
	if ($_=~/^>(.*)$/){
		my $id_old = $1;
#		print "$id_old\n";
		my @a=split /_/,$id_old;
		print ">$a[0]".'_'."$a[1]\n";
	}
	else{
		print "$_\n";
	}

}
