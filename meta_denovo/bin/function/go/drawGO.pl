#!/usr/bin/perl

#Author: yumingqi
#Date: 2014-11-01
#ediy by chenweitian 2016-12
use strict;
use Getopt::Long;
use FindBin '$Bin';
use File::Basename;
use File::Path 'mkpath';

my ($diff, $goprefix, $goclass, $outprefix, $Rscript, $convert);
GetOptions(
	"list=s" => \$diff,
	"goprefix=s" => \$goprefix,
	"goclass=s" => \$goclass,
	"outprefix:s" => \$outprefix,
	"Rscript:s" => \$Rscript,
	"convert:s" => \$convert
);
die <<USAGE unless $diff and $goprefix and $goclass and $outprefix;
Draw GO Classification of Differentially Expressed Genes of Sample.
	*-list        <file>    list of genes, format: <GeneID SomethingOrNot>
	*-goprefix    <path>    dir path which contains go files: *.C *.F *.P 
	*-goclass     <file>    go.class file
	 -outprefix   <prefix>  output prefix
	 -Rscript     <app.>    Rscript path
	 -convert     <app.>    convert path
2014-11-01 by Yu Mingqi
USAGE

$Rscript ||= "Rscript";
$convert ||= "convert";

my %geneID = map {(split)[0] => 1} `cat $diff`;
my (%gene2go, %go2content, %gene2content, %go2gene); # use goID as a bridge, connect gene to GO function
for ('C', 'F', 'P') {
	for (`cat $goprefix.$_`) {
		my ($id, $go) = (split)[1, 4];
		$gene2go{$id}{$go} = 1 if $geneID{$id}
	}
}
for (`cat $goclass`) {
	my ($a, $b, $go) = split /\t/;
	$go2content{$go}{"$a\t$b\n"} = 1
}

# open G2GO,">$outprefix.Gene2GO.xls" or die $!; print G2GO "Unigene\tGO ID\n";
# open GO2G,">$outprefix.GO2Gene.xls" or die $!; print GO2G "Ontoloty\tGO term\tNumber of Genes\tGene members\n";
open G2GO,">$outprefix.Gene2GO.txt" or die $!; print G2GO "基因\tGO编号\n";
open GO2G,">$outprefix.GO2Gene.txt" or die $!; print GO2G "本体(Ontoloty)\tGO term\t基因数\t基因\n";

for my $id (keys %gene2go) {
	my $out = $id;
	for my $go (keys %{$gene2go{$id}}) {
		$out .= "\t$go";
		for my $str (keys %{$go2content{$go}}) {
			$gene2content{$id}{$str} = 1;
			chomp($str); my @st = split /\t/,$str;
			$go2gene{$st[0]}{$st[1]}{$id} = 1;
# A gene may have 2 goID which both correspond to a same go class. So use hash reduce duplication.
		}
	}
	print G2GO "$out\n";
}

for my $cfp (sort keys %go2gene) {
	for my $cfp2 (sort keys %{$go2gene{$cfp}}) {
		my $num = 0;
		my $genes = "";
		for my $g (sort keys %{$go2gene{$cfp}{$cfp2}}) {
			$num++;
			$genes .= "$g;";
		}
		chop $genes;
		print GO2G "$cfp\t$cfp2\t$num\t$genes\n";
	}
}

close G2GO;
close GO2G;

die "Error: GeneIDs in -diff don't match GeneIDs in -goprefix. Please check.\n"	if (keys %gene2go) == 0;
die "Error: GOids in -goclass don't match GOids in -goprefix. Please check.\n"	if (keys %gene2content) == 0;

my %max_hash;
open MATRIX, ">$outprefix.matrix" or die $!;
for my $id (keys %gene2content) {
	for my $str (keys %{$gene2content{$id}}) {
		print MATRIX "$str\n";
		$max_hash{$str} ++;
	}
}
close MATRIX;
my @ys = sort {$b<=>$a} values %max_hash;
my $y_max = $ys[0] * 1.2;

open R, ">$outprefix.R" or die $!;
print R <<CODE;
pdf('$outprefix.GO.pdf', width=10, height=10)
library(ggplot2)
library(RColorBrewer)
library(grid)
d=read.table('$outprefix.matrix', sep='	')
df=data.frame(table(d\$V2, d\$V1)) # a new matrix
df\$Freq[df\$Freq==0]=NA # value = 0 -> value = NA
df=na.omit(df) # filter NA, for more beautiful pdf, but not necessary
df\$Var1=factor(df\$Var1, levels=rev(df\$Var1)) # make labels of GO func. as wish
colours=colorRampPalette(brewer.pal(9,"Set1"))(9)[1:3] 
order1<-df[order(-df\$Freq),]
order2<-order1[order(order1\$Var2),]
order2\$Var1<-factor(order2\$Var1,levels =order2\$Var1)
p<-ggplot(df,aes(x=Var1,y=Freq,fill=Var2)) + scale_fill_manual(values=colours)+
geom_bar(position="stack",stat="identity") +
theme(panel.background=element_rect(fill=NA,colour="grey"), panel.grid=element_line(color='grey'), panel.border=element_rect(fill='transparent',color='black'),axis.title=element_text(size=20), axis.text.x=element_text(color='black', size=10),axis.text.y=element_text(color='black', size=10))+
coord_flip()+labs(x='', y='Number of Genes')+theme(legend.title=element_blank(), legend.text=element_text(angle=270),legend.key.width=unit(2, 'mm'), legend.key.height=unit(5, 'cm'),legend.text.align=0.5,plot.title=element_text(face='bold', size=20))+scale_x_discrete(limits = rev(levels(order2\$Var1)))+scale_y_continuous(limits=c(0,$y_max), trans='sqrt')+geom_text(aes(label=Freq), hjust=-0.5, vjust=0.5, size = 3)
p
dev.off()
q()
CODE
close R;
`$Rscript $outprefix.R && $convert -density 300 -resize 30% $outprefix.GO.pdf $outprefix.GO.png`;
`rm $outprefix.matrix $outprefix.R`;
