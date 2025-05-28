#!/usr/bin/env perl 
# name: ttf2afm.pl
# func: get afm info from ttf file and saved in afm file
# require: PostScript::Font and Font::TTF
# if you want turn other format to afm, the web site:
# <http://www.zhuan-huan.com/font-converter.php>
# maybe useful.

use strict;
use warnings;
use utf8;
use PostScript::PrinterFontMetrics;
use PostScript::Resources;
use PostScript::Font::TTtoType42;
use File::Basename qw(basename);

die "Usage: perl $0 <path>\n" unless @ARGV == 1;

my $path = shift || '.';
my @files = _glob_files($path,'otf');

foreach my$file (@files)
{
	my $font = PostScript::Font::TTtoType42->new("$file");
	my $fname = basename($file);
	$fname = $1 if ($fname =~ /^(.+)\./);
	$font->write_afm("$fname.afm");
}

# get the files who has the specific suffix
sub _glob_files
{
	my $path = shift;
	my $suffix = shift;
	my @files;

	opendir DIR,$path or die "$path $!";
	while(my $filename = readdir(DIR))
	{
		push @files , "$path/$filename" if ($filename =~ /\.$suffix$/i)	
	}
	closedir DIR;
	return @files;
}
