#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin '$Bin';
use File::Basename qw(dirname basename);
use Cwd qw(abs_path);

sub usage{
        print STDERR"
        Draw COG Classification of Sample.
        *-catalog     <file>    *.class.catalog file
	 -outdir      <path>    outdir
	 -sample      <str>     sample name
         -Rscript     <app.>    Rscript path
         -convert     <app.>    convert path
         -help                  print this usage
        \n";
        exit;
}

my ($Catalog, $Outdir, $Rscript, $Convert, $Sam);
my $Help;
GetOptions(
        "catalog:s"  => \$Catalog,
        "outdir:s"   => \$Outdir,
        "sample:s"   => \$Sam,
        "Rscript:s"  => \$Rscript,
        "convert:s"  => \$Convert,
        "help"       => \$Help
);

&usage() if(!$Catalog || $Help);

$Rscript ||= "Rscript";
$Convert ||= "convert";
$Sam ||= "samplename";
$Outdir ||= '.';
mkdir $Outdir unless(-d $Outdir);
$Outdir = abs_path($Outdir);

my $prefix = "$Outdir/$Sam.COG";

my $y_max = 0;
open IN, "$Catalog" or die $!;
open MATRIX, "> $prefix.matrix" or die $!;
while(<IN>){
	chomp;
	next if(/^#/ || /^\s*$/);
	next if(/^Code/i);
	my @items = split /\t/;
	$y_max = ($y_max < $items[2])?$items[2]:$y_max;
	print MATRIX "$items[1]\t$items[2]\n";
}
close IN;
close MATRIX;
$y_max *= 1.2;

#my $color_numb = 24;
#my $color = 'colorRampPalette(brewer.pal(12,"Set3"))('.$color_numb.')';
### R ###
open R, ">$prefix.R" or die $!;
print R <<CODE;
pdf('$prefix.pdf', width=10, height=10)
d=read.table('$prefix.matrix', sep='\t')
#color = rainbow(24)
library(ggplot2)
library(grid)
#library(RColorBrewer)
NUM=length(d\$V2)
color=rep("#09BFFE",NUM)
ggplot(data=d,aes(x=V1,y=V2,fill=V1))+geom_bar(position='stack', stat='identity')+
     labs(title='COG Function Classification', x='', y='Number of Genes')+
     theme(panel.background = element_rect(fill='transparent'),
           panel.grid=element_line(color='grey'),
           panel.border=element_rect(fill='transparent',color='black'),
           legend.position="none", axis.text=element_text(color='black', size=12))+
     geom_text(aes(label=V2), hjust=-0.5, vjust=0.5, size = 4)+
     scale_y_continuous(limits=c(0,$y_max),trans='sqrt')+coord_flip()+
     scale_fill_manual(values=color)
dev.off()
CODE
close R;
`$Rscript $prefix.R && $Convert -density 300 -resize 30% $prefix.pdf $prefix.png`;


