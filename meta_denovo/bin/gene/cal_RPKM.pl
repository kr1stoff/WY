use strict;

my $list = shift;
my $count = shift;

my $total_reads;
open A,$list || die $!;
while(<A>){
	chomp;
	my @a = split /\t/,$_;
	$total_reads = $a[1];
}
close A;

print "Contig\tStart\tEnd\t基因\tReads_count\t基因长度\tTotal_reads\t基因丰度\n";

open B,$count || die $!;
while(<B>){
	chomp;
	my @b = split /\t/,$_;
	my $id = (split /_/,$b[0])[0];
	my $len = $b[2] - $b[1];
	my $rpkm = $b[4]*1000000000/($total_reads*$len);
	$rpkm = sprintf "%.2f",$rpkm;
	print "$_\t$len\t$total_reads\t$rpkm\n";
}
