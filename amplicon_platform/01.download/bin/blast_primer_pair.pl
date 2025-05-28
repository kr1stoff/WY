use strict;
use Data::Dumper;

my $table = shift;
my $in = shift;
my $name = shift;

my ($primer1, $primer2);
open T, $table || die $!;
while (<T>) {
    chomp;
    my ($n, $p1, $p2) = (split /\t/, $_)[0, 1, 2];
    if ($n eq "$name") {
        $primer1 = $p1;
        $primer2 = $p2;
    }
}
close T;

my (%hash_posi_positiveR, %hash_posi_positiveF, %hash_posi_negativeR, %hash_posi_negativeF);
open A, $in || die $!;
while (<A>) {
    chomp;
    next if /query/;
    my ($query, $seqid, $mismatch, $gap, $qstart, $sstart, $send, $slen, $taxid, $strand) = (split /\t/, $_)[0, 1, 4, 5, 6, 8, 9, 13, 15, 17];

    if ($strand eq "plus") {
        my $sstart_new;
        if ($sstart - $qstart + 1 <= 0) {
            $sstart_new = 1;
        }
        else {
            $sstart_new = $sstart - $qstart + 1; # 最前面上没比对上的延长，最后面没比对上的已经被过滤了
        }

        if ($query eq "$primer2") {
            # PRIMER-2
            $hash_posi_positiveR{$seqid} .= "$query:$sstart_new:$send:$mismatch:$gap:$qstart:$taxid" . "__";
        }
        else {
            # PRIMER-1
            $hash_posi_positiveF{$seqid} .= "$query:$sstart_new:$send:$mismatch:$gap:$qstart:$taxid" . "__";
        }
    }
    elsif ($strand eq "minus") {
        my $send_new;
        if ($sstart + $qstart - 1 > $slen) {
            $send_new = $slen;
        }
        else {
            $send_new = $sstart + $qstart - 1;
        }
        if ($query eq "$primer2") {
            # PRIMER-2
            my $sstart_new = $send;
            $hash_posi_negativeR{$seqid} .= "$query:$sstart_new:$send_new:$mismatch:$gap:$qstart:$taxid" . "__";
        }
        else {
            # PRIMER-1
            my $sstart_new = $send;
            $hash_posi_negativeF{$seqid} .= "$query:$sstart_new:$send_new:$mismatch:$gap:$qstart:$taxid" . "__";
        }
    }
}
close A;

print "#name\tid\ttaxid\tfr\tstrand1\tstart1\tend1\tmismatch1\tgap1\tstart_gap1\tfr\tstrand2\tstart2\tend2\tmismatch2\tgap2\tstart_gap2\tall_len\ttotal_len\tstrand\n";

foreach my $id (sort keys %hash_posi_positiveF) {
    if (exists $hash_posi_negativeR{$id}) {
        # print "$id\t$hash_posi_positiveF{$id}\t$hash_posi_negativeR{$id}\n";
        my @lists1 = (split /__/, $hash_posi_positiveF{$id});
        my @lists2 = (split /__/, $hash_posi_negativeR{$id});
        foreach my $posi1 (@lists1) {
            my ($primer1, $start1, $end1, $mismatch1, $gap1, $qstart1, $taxid) = (split /:/, $posi1)[0, 1, 2, 3, 4, 5, 6];
            foreach my $posi2 (@lists2) {
                # print "Test\t$posi1\t$posi2\n";
                my ($primer2, $start2, $end2, $mismatch2, $gap2, $qstart2) = (split /:/, $posi2)[0, 1, 2, 3, 4, 5];
                # print "Posi\t$start2\t $start1\n";
                if ($start2 - $start1 < 2000 && $start2 - $start1 > 0) {
                    my $len1 = $end1 - $start1;           # F primer length
                    my $len2 = $end2 - $start2;           # R primer length
                    my $len_insert = $start2 - $end1 + 1; # insert length
                    my $total_len = $len1 + $len_insert + $len2;
                    my $start_gap1 = $qstart1 - 1;
                    my $start_gap2 = $qstart2 - 1;
                    # print "$id\t$primer1\t+\t$start1\t$end1\t$mismatch1\t$gap1\t$primer2\t-\t$start2\t$end2\t$mismatch2\t$gap2\n";
                    print "$name\t$id\t$taxid\tF\t+\t$start1\t$end1\t$mismatch1\t$gap1\t$start_gap1\tR\t-\t$start2\t$end2\t$mismatch2\t$gap2\t$start_gap2\t$len1+$len_insert+$len2\t$total_len\t+\n";
                }
            }
        }
    }
}

# print "#############################\n";

foreach my $id (sort keys %hash_posi_positiveR) {
    if (exists $hash_posi_negativeF{$id}) {
        # print "$id\t$hash_posi_positiveR{$id}\t$hash_posi_negativeF{$id}\n";
        my @lists1 = (split /__/, $hash_posi_positiveR{$id});
        my @lists2 = (split /__/, $hash_posi_negativeF{$id});
        foreach my $posi1 (@lists1) {
            my ($primer1, $start1, $end1, $mismatch1, $gap1, $qstart1, $taxid) = (split /:/, $posi1)[0, 1, 2, 3, 4, 5, 6];
            foreach my $posi2 (@lists2) {
                # print "Test\t$posi1\t$posi2\n";
                my ($primer2, $start2, $end2, $mismatch2, $gap2, $qstart2) = (split /:/, $posi2)[0, 1, 2, 3, 4, 5];
                # print "Posi\t$start2\t $start1\n";
                if ($start2 - $start1 < 2000 && $start2 - $start1 > 0) {
                    my $len1 = $end1 - $start1;           # R primer length
                    my $len2 = $end2 - $start2;           # F primer length
                    my $len_insert = $start2 - $end1 + 1; # insert length
                    my $total_len = $len1 + $len_insert + $len2;
                    my $start_gap1 = $qstart1 - 1;
                    my $start_gap2 = $qstart2 - 1;
                    # print "$id\t$primer1\t+\t$start1\t$end1\t$mismatch1\t$gap1\t$primer2\t-\t$start2\t$end2\t$mismatch2\t$gap2\n";
                    print "$name\t$id\t$taxid\tR\t+\t$start1\t$end1\t$mismatch1\t$gap1\t$start_gap1\tF\t-\t$start2\t$end2\t$mismatch2\t$gap2\t$start_gap2\t$len1+$len_insert+$len2\t$total_len\t-\n";
                }
            }
        }
    }
}


