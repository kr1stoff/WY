use strict;

my $chromosome = shift;
my $in = shift;
my $out_chromosome = shift;
my $out_in = shift;
my $out_stat = shift;

my %hash;
my @all_gcf;
my $taxid;
open T,$chromosome || die $!;
while(<T>){
    chomp;
    next if /^taxon_id/;
    $taxid = (split /\t/,$_)[0];
    my ($name,$all_id) = (split /\t/,$_)[1,2];
    push @all_gcf,$name;
    my @ids = split /,/,$all_id;
    foreach my $id (@ids){
        $hash{$id} = $name;
    }
}
close T;

my @chrs;
my $primer_id;
open A,$in || die $!;
open OUT1,">$out_in" || die $!;
while(<A>){
    chomp;
    next if /^#/;
    $primer_id = (split /\t/,$_)[0];
    my $name = (split /\t/,$_)[1];
    my $chr = $hash{$name};
    print OUT1 "$_\t$chr\n";
    push @chrs,$chr;
}
close A;
close OUT1;

my %chr;
foreach (@chrs) {$chr{$_}++;}

open OUT2,">$out_chromosome" || die $!;

my $all_gcf_num = @all_gcf;
my $epcr_num = 0;
foreach my $gcf (@all_gcf){
    if (exists $chr{$gcf}){
        $epcr_num  += 1;
    }
    print OUT2 "$taxid\t$gcf\t$chr{$gcf}\n";
}

my $per = sprintf "%.2f",(100*$epcr_num/$all_gcf_num);

open OUT_stat,">$out_stat" || die $!;
print OUT_stat "#Name\tPercent\tMatch_num\tAll_num\n";
print OUT_stat "$primer_id\t$per\t$epcr_num\t$all_gcf_num\n";
close OUT_stat;