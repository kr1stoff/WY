use strict;

my $table = shift;
my $blast = shift;
my $outdir = shift;

$outdir ||= "./";

my %hash;
open B,$blast || die $!;
while(<B>){
    chomp;
    my $id = (split /\t/,$_)[0];
    $hash{$id} .= $_."\n";
}
close B;

open A,$table || die $!;
while(<A>){
    chomp;
    my ($id,$f,$r) = (split /\t/,$_)[0,1,2];
    `mkdir -p $outdir/$id`;
    open OUT,">$outdir/$id/$id.blast" || die $!;
    print OUT "$hash{$f}";
    print OUT "$hash{$r}";
    open OUT2,">$outdir/$id/$id.input.table" || die $!;
    print OUT2 "$id\t$f\t$r\n";
}
close A;
