use strict;

my $lineage = shift;
my $abun = shift;

my %hash;
open A,$lineage || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    $hash{$a[0]} = $a[1];
}
close A;

print "Contig ID\tContig丰度\tContig分类\n";
open B,$abun || die $!;
while(<B>){
    chomp;
    my @b = split /\t/,$_;
    my $abun = sprintf "%.2f",$b[1];
    if (exists $hash{$b[0]}){
        print "$b[0]\t$abun\t$hash{$b[0]}\n";
    }
}
close B;