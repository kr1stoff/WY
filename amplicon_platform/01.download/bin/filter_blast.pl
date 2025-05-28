use strict;

my $blast = shift;

# print "query\tseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovs\ttaxid\tsscinames\tsstrand\n";

open A, $blast || die $!;
while (<A>) {
    chomp;
    my ($qaccver, $saccver, $mismatch, $gapopen, $qstart, $qend, $qlen) = (split /\t/, $_)[0, 1, 4, 5, 6, 7, 12];
    my $all_mis = $qlen - ($qend - $qstart + 1) + $mismatch + $gapopen;
    if ($qend == $qlen && $all_mis <= 2) {
        print "$_\n";
    }
}