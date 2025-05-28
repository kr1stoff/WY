use strict;

my $in = shift;

my (%gene,%pathway,%K);
my $total;

open A,$in || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    $pathway{$a[0]} = "$a[1]\t$a[2]\t$a[3]";
    $gene{$a[0]} .= $a[5].";";

    my @tmp = split /;/,$a[5];
    $total += ($#tmp+1);
    $K{$a[0]} .= $a[4].";";
}
close A;

print "KEGG_A_class\tKEGG_B_class\tPathway\tCount ($total)\tPathway ID\tGenes\tK_IDs\n";

foreach my $ko (sort keys %gene){
    my @all = split /;/,$gene{$ko};
    print "$ko\t$pathway{$ko}\t$#all\t$gene{$ko}\t$K{$ko}\n";
}
