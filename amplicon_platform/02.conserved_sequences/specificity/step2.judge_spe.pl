use strict;

my $list = shift;
my $blast6 = shift;

my %hash;
open A,$list || die $!;
while(<A>){
    chomp;
    $hash{$_} = 1;
}
close A;

open B,$blast6 || die $!;
while(<B>){
    chomp;
    my $taxid = (split /\t/,$_)[17];
    if (exists $hash{$taxid}){
        print "$_\tSame\n";
    }
    else{
        print "$_\tDiff\n";
    }
}