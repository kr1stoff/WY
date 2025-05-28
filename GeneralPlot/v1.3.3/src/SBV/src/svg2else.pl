#!/usr/bin/perl -w
# name: svg2else.pl
# func: turn svg format file to other format by 'inkscape'

use strict;
use warnings;

use Getopt::Std;

my %opts = (f=>'png',A=>'draw',d=>300,p=>"D:/Program Files/Inkscape/inkscape");
getopts('f:d:A:p:h',\%opts);

&usage() if (@ARGV == 0 || $opts{h});

my $cmd = "";
# check the inkscape executable path
my $inkscape = $opts{p};
#die "the inkscape executable path is not exists,\n[$inkscape]" if (! -e $inkscape);
$cmd .= "\"$inkscape\" ";

# the output file format, default is png
# now inkscape suppor format: png(-e), 
my %format = (png=>'-e',ps=>'-P',eps=>'-E',pdf=>'-A',emf=>'-M');
$opts{f} = lc $opts{f};
die "the format is not supprted in inkscape, [$opts{f}]" if (! exists $format{$opts{f}});

# set the dpi, default is 300
$cmd .= "-d $opts{d} ";

# the output area
# draw: exported area is the entire drawing 
# page: exported area is the entire page
# default is draw
my %area = (draw=>'-D',page=>'-C');
if (! exists $area{$opts{A}})
{
	warn ("the output area is either 'draw' or 'page', [$opts{A}]");
	$opts{A} = 'draw';
}
$cmd .= $area{$opts{A}} . " ";

# set the export background color as white
$cmd .= "-b white ";

# the svg files need to be turn 
my @files = @ARGV;

foreach my$file (@files)
{
	die "this file is not in svg format, [$file]" if ($file !~ /\.svg$/i);
	my $fname = $1 if ($file =~ /(.+)\.svg/i);
	my $outname = "$fname.$opts{f}";
	my $newcmd = "$cmd $format{$opts{f}} $outname $file";
	system($newcmd);
}


sub usage
{
	print <<HELP;
Usage:   perl $0 [options] <*.svg>[s]

Options: -f STR    the export out file format,png|pdf|emf|ps|eps [png]
         -d INT    the dpi of out file, [300]
         -A STR    the export area of the svg file, draw|page [draw]
         -p STR    the path of the inkscape, you can set it in line 10
         -h        show this information

Note:    the default output background color is white in line 43, 
         for more function, you should use the inkscape directly.
HELP
	exit;
}