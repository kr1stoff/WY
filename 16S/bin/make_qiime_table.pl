#!/usr/bin/perl -w
use strict;
die "perl $0 <all.otus.shared|otus.tb.tsv> <tax_assignments.txt>\n" unless (@ARGV ge 1);

###### read otu table ######
my (%otu2tags, $tag, @groups, @otus);

if ($ARGV[0] =~ /\.shared$/) {
    open IN, "<$ARGV[0]" or die "Can't open the file! $!";
    my $head = <IN>;
    chomp $head;
    (undef, undef, undef, @otus) = split "\t", $head;
    while (<IN>) {
        chomp;
        my (undef, $group, undef, @tags) = split "\t", $_;
        push @groups, $group;
        for my $i (0 .. $#otus) {
            push @{$otu2tags{$otus[$i]}}, $tags[$i];
        }
    }
    close IN;
}
else {
    open IN, "<$ARGV[0]" or die "Can't open the file! $!";
    my $head = <IN>;
    chomp $head;
    ($tag, @groups) = split "\t", $head;
    while (<IN>) {
        chomp;
        my ($otu, @tags) = split "\t", $_;
        push @otus, $otu;
        @{$otu2tags{$otu}} = @tags;
    }
}

###### read align table ######
if ($ARGV[1] =~ /\.report$/) {
    my %otu2hit;
    open IN, "<$ARGV[1]" or die "Can't open the file! $!";
    <IN>;
    while (<IN>) {
        my ($otu, undef, $hit, undef) = split "\t", $_, 4;
        $otu2hit{$otu} = $hit;
    }
    print "OTU ID\t" . join("\t", @groups) . "\n";
    my %hit2tags;
    for my $otu (@otus) {
        my $hit = $otu2hit{$otu};
        if (exists $hit2tags{$hit}) {
            for my $i (0 .. $#groups) {
                ${$hit2tags{$hit}}[$i] += ${$otu2tags{$otu}}[$i];
            }
        }
        else {
            @{$hit2tags{$hit}} = @{$otu2tags{$otu}};
        }
    }
    for my $hit (keys %hit2tags) {
        print "$hit\t" . join("\t", @{$hit2tags{$hit}}) . "\n";
    }
    exit 0;
}

open IN, "<$ARGV[1]" or die "Can't open the file! $!";
<IN>;
my (%otu2tax, %otu2tax2);
while (<IN>) {
    chomp;
    my ($otu, $tax, undef) = split "\t", $_;
    ## tax for annot
    $tax = "Root" if $tax =~ /No blast hit/;
    $tax =~ s/Root;//;
    $tax =~ s/\ /_/g;
    $tax =~ s/o__Chloroplast.*/o__Chloroplast/; # add 20210111
    $tax =~ s/__\[/__/g;
    $tax =~ s/\];/;/g;
    $otu2tax2{$otu} = $tax;

    ## tax for function
    $tax =~ s/;[a-z]__/;/g;
    $tax =~ s/_sp_/_sp._/;
    $tax =~ s/_/\ /g;
    $tax =~ s/;*$/;/;
    $tax =~ s/_sp;/_sp.;/;
    $tax = ";" if $tax eq "";
    $otu2tax{$otu} = $tax;
}
close IN;

## make qiime otu table ##
my $label = $tag;
print "$label\t", join("\t", @groups), "\tTaxonomy\n";
for my $i (0 .. $#otus) {
    my $otu = $otus[$i];
    print "$otu\t", join("\t", @{$otu2tags{$otu}}), "\t$otu2tax2{$otu}\n";
}

