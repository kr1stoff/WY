use strict;
use Data::Dumper;

my $step3 = shift;
my $stat = shift;

my %hash;
open A,$step3 || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    $hash{$a[0]} = "$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]";
}
close A;

print "Name\tCon_count\tCon_total\tCon_percent\tCon_min_iden\tCon_min_cov\tCon_av_iden\tCon_av_cov\tSpe_per_same\tSpe_diff_max\tSpe_diff_min\tSpe_same_max\tSpe_same_min\tSpe_total\n";
open C,$stat || die $!;
while(<C>){
    chomp;
    next if /query_id/;
    my @t = split /\t/,$_;
    my $name = $t[0];
    # print "$name\n";
    if (exists $hash{$name}){
        print "$_\t$hash{$name}\n";
    }
}
close C;