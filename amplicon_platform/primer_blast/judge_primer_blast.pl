use strict;

my $list = shift;
my $blast6 = shift;
my $line = shift;

$line ||= "16";

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
    next if /^#/;
    my $taxid = (split /\t/,$_)[$line-1];
    if (exists $hash{$taxid}){
        print "$_\tSame\n";
    }
    else{
        print "$_\tDiff\n";
    }
}