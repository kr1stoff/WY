use strict;
use List::Util qw(max min sum);
use Data::Dumper;

my $blast = shift;

my (%hash_count,%hash_same,%hash_diff);
open A,$blast || die $!;
while(<A>){
	chomp;
	my @a = split /\t/,$_;
	my ($gene,$identity,$compare) = @a[0,2,20];
    $hash_count{$gene} += 1;
    if ($compare eq "Same"){
        $hash_same{$gene} .= $identity."__";
    }
	if ($compare eq "Diff"){
        $hash_diff{$gene} .= $identity."__";
    }
}
# print Dumper (\%hash_count);

print "Gene\tPer_same\tDiff_max\tDiff_min\tSame_max\tSame_min\tTotal\n";
foreach my $gene (sort keys %hash_count){
    my $total = $hash_count{$gene};
    if (exists $hash_same{$gene}){
        # next if ! exists $hash_same{$gene};
        my @b = split /__/,$hash_same{$gene};
        my $num_same = @b;
        my $per_same = $num_same*100/$total;
        $per_same=sprintf "%.2f",$per_same;

        my ($diff_max,$diff_min);
        my ($same_min,$same_max);

        $same_min=min @b;
        $same_max=max @b;
        if (exists $hash_diff{$gene}){
            my @d = split /__/,$hash_diff{$gene};
            $diff_min=min @d;
            $diff_max=max @d;
        }
        else{
            $diff_max = 85;
            $diff_min = 85;
        }
        print "$gene\t$per_same\t$diff_max\t$diff_min\t$same_max\t$same_min\t$total\n";
    }
    else{
        print "$gene\tAll_fail\n";
    }

}
close A;