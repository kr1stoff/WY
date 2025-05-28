#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use FindBin qw($Bin $RealBin);
use YAML qw(LoadFile);

# 读取所用软件
my $f_software = "$RealBin/../config/software.yml";
our $app = LoadFile($f_software);
my $perl = $app->{"perl"};

my $example1 = "perl $0 DataFrame 1,2,3 /Bio/Bin/pipeline/GeneralPlot/v1.0/bin/general_plot.pl bar -outprefix output -barfill T";
my $example2 = "perl $0 Tab2Transpose t /Bio/Bin/pipeline/GeneralPlot/v1.0/bin/general_plot.pl bar -outprefix output -barfill T";
my $message = "perl $0 <file> <1,2,3 or t> <general_plot> args\nexample:\n$example1\n$example2\n";
die $message unless (@ARGV > 3);

my @cols = split ',', $ARGV[1];
my @samples;
my @variables;
my %hash;
my %hash2;
open IN, "<$ARGV[0]" or die "Can't open the file! $!";
my $head = <IN>;
if ($ARGV[1] =~ /,/) {
    die "three columns choose need!" if @cols != 3;
    while (<IN>) {
        chomp;
        my @tmps = split "\t", $_;
        my $v1 = $tmps[$cols[0] - 1];
        my $v2 = $tmps[$cols[1] - 1];
        my $v3 = $tmps[$cols[2] - 1];
        push @samples, $v1 unless exists $hash{$v1};
        push @variables, $v2 unless exists $hash2{$v2};
        $hash{$v1}{$v2} = $v3;
        $hash2{$v2}++;
    }
}
elsif ($ARGV[1] eq 't') {
    chomp $head;
    (undef, @samples) = split "\t", $head;
    while (<IN>) {
        chomp;
        my ($id, @tmps) = split "\t", $_;
        push @variables, $id;
        map {$hash{$samples[$_]}{$id} = $tmps[$_]} 0 .. $#tmps;
    }
    @variables = reverse(@variables);
}
else {
    die "[Error] Unknown parameter ($ARGV[1])!\n";
}
close IN;

my $data4fig = "$ARGV[0].forFig.tmp";
open OUT, ">$data4fig" or die "Can't write the file! $!";
print OUT "var\t", join("\t", @variables), "\n";
for my $sam (@samples) {
    print OUT $sam;
    for my $var (@variables) {
        $hash{$sam}{$var} ||= 0;
        print OUT "\t", $hash{$sam}{$var};
    }
    print OUT "\n";
}
close OUT;

my $cmd = "$perl " . join(" ", @ARGV[2 .. $#ARGV]) . " -file $data4fig";
print "[CMD] $cmd\n";
my $sig = system($cmd);
print("[LOG] Done with sig $sig\n");
system("rm -rf $data4fig") if $sig == 0;

