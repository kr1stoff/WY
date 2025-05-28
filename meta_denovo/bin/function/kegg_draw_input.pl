use strict;
use Data::Dumper;

my $kopath = shift;
my $format = shift;

my %hash;
open A,$kopath || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    my $K = (split /,/,$a[1])[0];
    $hash{$K} .= $a[0].";";
}
close A;
# print Dumper (\%hash);

#KEGG_A_class    KEGG_B_class    Pathway Count (2969)    Pathway ID      Genes   K_IDs

open B,$format || die $!;
while(<B>){
    chomp;
    my @b = split /\t/,$_;
    if (exists $hash{$b[4]}){
        print "$_\t$hash{$b[4]}\n";    
    }
}
close B;
