#!/usr/bin/perl
#########################################################
# Author: handsye
# Created Time : Mon 13 Jun 2022 04:07:29 PM CST
# File Name: result.pl
# Version: v0.1.0
#########################################################

use strict;
use warnings;
use Data::Dumper;

my @ids;
my %hash;
open IN, "<$ARGV[0]" or die $!;
while (<IN>) {
    chomp;
    next if $. == 1;
    my @line = split /\s+/, $_;
    push @ids, $line[0];
}
close IN;

open IY, "<$ARGV[0]" or die $!;
while (<IY>) {
    chomp;
    next if $. == 1;
    my @line = split /\s+/, $_;
    if ($line[0] eq $ARGV[1]) {
        for (my $i = 1; $i < scalar(@ids); $i = $i + 1) {
            $hash{$line[$i]} = $ids[$i - 1];
        }
    }
    else {
        next;
    }
}
close IY;

my @key = (sort keys %hash);
my $p = $hash{$key[1]};
$p =~ s/\.\w//g;
my $type;
open TY, "<$ARGV[2]" or die $!;
while (<TY>) {
    chomp;
    my @line = split /\t/, $_;
    if ($line[1] =~ $p) {
        $type = $line[0];
    }
}
close TY;
open OU, ">$ARGV[3]/result.tsv" or die $!;
print OU "ID\tTyping\tNeighbor\n";
print OU "$ARGV[1]\t$type\t$hash{$key[1]}\n";

