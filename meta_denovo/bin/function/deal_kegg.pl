use strict;

my $kotxt = shift;


open A,$kotxt || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    if ($a[5] =~ /(.*)\[PATH:(.*)\]/){
        my $pathway = $1;
        my $ko = $2;
        print "$ko\t$a[1]\t$a[3]\t$pathway\t$a[6]\n";
    }
}
close A;
