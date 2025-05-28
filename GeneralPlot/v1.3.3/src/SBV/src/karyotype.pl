#! /usr/bin/perl -w
#name: karyotype.pl

use strict;

die "Usage: perl $0 <file> [color]\n" 
   unless (@ARGV==1 || @ARGV==2);

my $file = shift;
my $color = shift || "grey";

len2KARYO($file,$color);

# input: chr length list file (2 fileds: chr name<TAB>length)
# return: karyotype
sub len2KARYO
{
    my $file = shift;
    my $color = shift;
    
    open FH,$file || die "$file $!";
    while(<FH>)
    {
        next if (/^#/);
        chomp;
        my ($chrName,$len) = split;
        
        print "chr\t-\t$chrName\t$chrName\t0\t$len\t$color\n";
    }
    close FH;
}
