#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: drawSymbols.pl
#
#        USAGE: ./drawSymbols.pl  
#
#  DESCRIPTION: draw All Symbols in SBV
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Peng Ai (), apanhui@gmail.com
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 09/17/2013 11:23:58 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/..";
use lib "$FindBin::RealBin/../lib";

use SVG;
use SBV::STONE::SYMBOL;

# init some parameter
my $margin = 20;
my $ox = $margin;
my $oy = $margin;
my $num = 44;
my ($len1,$len2,$len3) = (60,200,100);
my $x1 = $ox + $len1;
my $x2 = $x1 + $len2;
my $x3 = $x2 + $len3;
my $unitH = 30;
my $width = $len1 + $len2 + $len3;
my $height = $unitH * ($num + 1 - 2); # 38,39 are under-development

my $svg = SVG->new(width=>$width+2*$margin,height=>$height+2*$margin,id=>'symbols_legend_map');

# define the css
my $defs = $svg->defs();
my $css = <<CSS;
text {
	fill:#000;
	font-family:arial;
	font-size:12px;
}

text.header
{
	font-weight:bold;	
}

line 
{
	stroke-width:1;
	stroke:#000;
}

line.boundary 
{
	stroke:#fff;
}
CSS
$defs->style(type=>"text/css")->CDATA($css);

# set the background and column split line
$svg->rect(x=>$ox,y=>$oy,width=>$width,height=>$height,style=>"fill:#ccc;stroke-width:0");
$svg->line(x1=>$x1,x2=>$x1,y1=>$oy,y2=>$oy+$height,class=>"boundary");
$svg->line(x1=>$x2,x2=>$x2,y1=>$oy,y2=>$oy+$height,class=>"boundary");

# init font
my $font = SBV::Font->new("family:arial;size:12px;weight:normal;style:normal");
my $textH = $font->fetch_text_height;

# add the header info
# ID		shape
my $textY = $oy+$unitH/2+$textH/2;
$svg->text(x=>$ox+10,y=>$textY,class=>"header")->cdata("Code");
$svg->text(x=>$x1+10,y=>$textY,class=>"header")->cdata("Shape");
$svg->text(x=>$x2+10,y=>$textY,class=>"header")->cdata("Example");
$oy += $unitH;
$svg->line(x1=>$ox,x2=>$x3,y1=>$oy,y2=>$oy,class=>"boundary");

# init the shape desc
my @shape = ('rectangle(rect)','circle','diamond', # 0,1,2
'up pointing regular triangle','low pointing regular triangle', # 3,4
'cross','error','cross and error','rect and cross','rect and error', # 5,6,7,8,9
'email','circle and cross','circle and error','diamond and cross', # 10,11,12,13
'3 and 4','five points star','six points star','horizontal line', # 14,15,16,17
'horizontal line and a point','boxplot','vertical line','vertical line and a point',# 18,19,20,21
'box bar','err bar','left pointing pencil','right pointing pencil',# 22,23,24,25
'right pointing triangle','left pointing triangle', # 26,27
'right pointing pentagram','left pointing pentagram', # 28,29
'up pointing pentagram','down pointing pentagram', # 30,31
'horizontal hexagon','vertical hexagon', # 32,33
'octagon','1/3 height rect','ellipse','rounded rect', # 34,35,36,37
'conventional heart','implicit heart','left pointing arrow','right pointing arrow', # 38,39,40,41
'fan (Center on the left)','semicircle', # 42,43
);

my @colors = SBV::Colors::rainbow(44);
$colors[35] = "#000";

for my$i( 0 .. $num-1 )
{
	next if ($i == 38 || $i==39);
	$textY = $oy+$unitH/2+$textH/2;
	$svg->text(x=>$ox+10,y=>$textY)->cdata($i);
	$svg->text(x=>$x1+10,y=>$textY)->cdata($shape[$i]);
	
	my $sw = $i <= 23 ? $unitH-10     : 
             $i == 42 ? ($unitH-10)*2 : 
             $i == 43 ? ($unitH-10)*2 : $len3-20;
	my $usex = $x2 + $len3/2 - $sw/2;

	my $style;
	if (($i<=13 && $i >= 5) || $i == 23 )
	{
		$style = {color=>$colors[$i]};
	}
	else
	{
		$style = {fill=>$colors[$i]};
	}

	my $id = SBV::STONE::SYMBOL::new($i,width=>$sw,height=>$unitH-10,parent=>$defs,%$style);
	$svg->use(x=>$usex,y=>$oy+5,width=>$sw,height=>$unitH-10,-href=>"#$id");
	
	$oy += $unitH;
	$svg->line(x1=>$ox,x2=>$x3,y1=>$oy,y2=>$oy,class=>"boundary");
}

open OUT,">","symbols.svg";
print OUT $svg->xmlify;
close OUT;
