use strict;

my $name = shift;
my $input_table = shift;
my $blast_stat = shift;
my $pair_stat = shift;

open A, $input_table || die $!;
my ($f_name, $r_name);
while (<A>) {
    chomp;
    my ($n, $f, $r) = (split /\t/, $_)[0, 1, 2];
    if ($n eq $name) {
        $f_name = $f;
        $r_name = $r;
    }
}
close A;

open B, $blast_stat || die $!;
my ($f_per, $f_tol, $r_per, $r_tol);
while (<B>) {
    chomp;
    next if /^#/;
    my ($na, $p, $t) = (split /\t/, $_)[0, 1, 2];
    if ($na eq $f_name) {
        $f_per = $p;
        $f_tol = $t;
    }
    elsif ($na eq $r_name) {
        $r_per = $p;
        $r_tol = $t;
    }
}
close B;

open C, $pair_stat || die $!;
my $pair_per ||= "0";
my $pari_tol ||= '-';
while (<C>) {
    chomp;
    next if /^#/;
    my ($na, $p, $t) = (split /\t/, $_)[0, 1, 2];
    $pair_per = $p;
    $pari_tol = $t;
}
close C;

print "#Name\tPair_percent\tPair_total\tF_percent\tF_total\tR_percent\tR_total\n";
print "$name\t$pair_per\t$pari_tol\t$f_per\t$f_tol\t$r_per\t$r_tol\n";