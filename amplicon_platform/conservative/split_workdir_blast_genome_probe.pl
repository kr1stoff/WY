use strict;

my $table = shift;
my $fa = shift;
my $outdir = shift;

$outdir ||= "./";

my %hash_fa;
read_fasta($fa,\%hash_fa);

open A,$table || die $!;
while(<A>){
    chomp;
    my ($id,$f,$r) = (split /\t/,$_)[0,1,2];
    `mkdir -p $outdir/$id`;
    open OUT_FA,">$outdir/$id/$id.fa" || die $!;
    open OUT_TABLE,">$outdir/$id/$id.input.table" || die $!;
    print OUT_FA ">$f\n$hash_fa{$f}";
    print OUT_TABLE "$id\t$f\n";
    close OUT_FA;
    close OUT_TABLE;
}
close A;

sub read_fasta{
        my ($file,$bait_hp) = @_;
        open(IN,$file)||die("fail to open $file\n");
        $/=">";<IN>;$/="\n";
        while (<IN>) {
                my $title=$_;
                chomp $title;
                $/=">";
                my $seq=<IN>;
                chomp $seq;
                $/="\n";
                my @temp=split(/\s+/,$title);
                my $id = $temp[0];

                $$bait_hp{$id}=$seq;

        }
        close(IN);
}

