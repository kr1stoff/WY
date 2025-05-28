use strict;

my $all_primer = shift;
my $all_probe = shift;
my $table = shift;

my (%primer,%probe);
read_file2hash($all_primer,\%primer);
read_file2hash($all_probe,\%probe);

print "#Name\tPrimer_perecnt\tPrimer_total\tProbe_percent\tProbe_total\tGenome_num\n";
open A,$table || die $!;
while(<A>){
    chomp;
    my ($id,$taxid) = (split /\t/,$_)[0,-1];
    my $primer_out;
    if (exists $primer{$id}){
        my ($per,$num,$total_num) = (split /\t/,$primer{$id})[0,1,2];
        $primer_out = "$per\t$num";
    }
    else{
        $primer_out = "-\t-";
    }

    my $probe_out;
    if (exists $probe{$id}){
        $probe_out = $probe{$id};
    }
    else{
        $probe_out = "-\t-\t-";
    }

    print "$id\t$primer_out\t$probe_out\n";
}
close A;

sub read_file2hash{
    my ($file,$hash) = @_;
    open T,$file || die $!;
    while(<T>){
        chomp;
        my @line = split /\t/,$_;
        my $name = $line[0];
        $$hash{$name}  = join "\t",@line[1..$#line];
    }
    close T;
}

