use strict;

my $desc_table = shift;
my $in = shift;
my $class_desc = shift;
my $out_number = shift;

my %hash;
open A,$desc_table || die $!;
while(<A>){
	chomp;
	my @t = split /\t/,$_;
	$hash{$t[0]} = $t[1];
}

my %hash_num;
open B,$in || die $!;
open OUT1,">$class_desc" || die $!;
while(<B>){
	chomp;
	next if /^#/;
	my @a = split/\t/,$_;
	my @b = split //,$a[1];
	my @desc;
	for my $i (@b){
		$hash_num{$i}+=1;
		my $desc = $hash{$i};
		push @desc,$desc;
	}
	print OUT1 "$a[0]\t";
	print OUT1 join (";",@desc)."\n";
}
close B;
close OUT1;

open OUT2,">$out_number" || die $!;
print OUT2 "Code\tFunctional-Categories\tGene-Number\n";
foreach my $class (sort keys %hash_num){
	print OUT2 "$class\t$hash{$class}\t$hash_num{$class}\n";
}
close OUT2;
