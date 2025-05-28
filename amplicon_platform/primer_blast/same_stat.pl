use strict;
use List::Util qw(max min sum);
use Data::Dumper;

my $blast = shift;

my (%hash_count,%hash_same,%hash_diff);
open A,$blast || die $!;
while(<A>){
	chomp;
    next if /query/;
    next if /^#/;
	my @a = split /\t/,$_;
	my ($name,$compare) = @a[0,-1];
    $hash_count{$name} += 1;
    if ($compare eq "Same"){
        $hash_same{$name} += 1;
    }
}

print "#Name\tPer_same\tTotal\n";
foreach my $name (sort keys %hash_count){
    my $total = $hash_count{$name};
    if (exists $hash_same{$name}){
        # next if ! exists $hash_same{$name};
        my $num_same = $hash_same{$name};
        my $per_same = $num_same*100/$total;
        $per_same=sprintf "%.2f",$per_same;

        print "$name\t$per_same\t$total\n";
    }
    else{
        print "$name\t0\t$total\n";
    }
}
close A;