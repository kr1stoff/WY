#!usr/bin/perl -w
use strict;
die "perl $0 <*Lefse.res> <outdir>\n" if @ARGV == 0;

my ($res, $outdir) = @ARGV;
$outdir ||= ".";

open my $fh_res, "<", $res or die $!;
chomp(my $outfile = `basename $res`);
$outfile =~ s/\.res/\.Diff\.xls/;
open my $ofh, ">", "$outdir/$outfile" or die $!;

print {$ofh} "feature\tlog_hm\tclass_hm\tLDA score\tp value\n";

while (<$fh_res>) {
    chomp;
    my ($feature, $log_hm, $class_hm, $score, $pvalue) = split /\t/, $_;

    if ($pvalue ne "-" && $log_hm ne "" && $class_hm ne "" && $score ne "") {
        print $ofh "$_\n";
    }
}
close $fh_res;




