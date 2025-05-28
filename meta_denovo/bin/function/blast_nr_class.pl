#!/usr/bin/perl -w

=head1 Name

	blast_nr_class.pl  --  statistics of blastx nr result

=head1 Description
  This program is designed for doing statistics of blastx nr result.

=head1 Usage

  perl blast_nr_class.pl <options>

=head1 sample 
	perl blast_nr_class.pl -nr xxxx.fa.blast.Nr.xls -nohead -out ./
=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename;
use lib $Bin;
my ($nr, $nohead, $outdir, $help);
GetOptions(
	"nr:s"=>\$nr,
	"nohead"=>\$nohead,
	"outdir:s"=>\$outdir,
	"help!"=>\$help
);

die `pod2text $0` if($help);
die `pod2text $0` unless ($nr and $outdir);

# remove duplication
if(defined $nohead){
	system "sort -k 1,1 -u $nr > $nr.sort";
}
else{
	system "sed 1d $nr | sort -k 1,1 -u > $nr.sort";
}

our $total_gene = 0;
my ($gene_name, $nr_id, $identity, $evalue, $annotation);
my %ehash;	# evlue hash
my %ihash;	# identity hash
our $min_iden = 1;
my %shash;	# species hash

# read sorted nr blast result
open NR,"$nr.sort" or die $!;
print STDERR "Statistics nr information....\n";
while(<NR>){
        chomp;
        $total_gene += 1;
        my @nr_arr = split(/\t/, $_,13);
        ($gene_name, $nr_id, $identity, $evalue, $annotation) = ($nr_arr[0], $nr_arr[1], $nr_arr[2], $nr_arr[10], $nr_arr[12]);	
	species_func($annotation);
}
close NR;

# output the species info
our @species_array;
my $species_count = 0;
open SOUT,">$outdir/species_statistic.xls";
print SOUT "species\tgene numbers\tpercentage\n";
foreach my $sKey(sort {$shash{$b} <=> $shash{$a}} keys %shash){
        my $sPer = sprintf("%.2f", 100 * $shash{$sKey} / $total_gene);
	if($sKey ne "other"){
		printf SOUT "%s\t%d\t%.2f%s\n", $sKey, $shash{$sKey}, $sPer, "%";
		if($species_count < 4 && $sPer >= 1){
			$species_count += 1;
	                push @species_array, $sKey;
        	}
	        else{
        	        $shash{'other'} += $shash{$sKey};
		}
        }
}
#print STDERR $shash{'other'};
push @species_array, 'other';

# drow picture
draw_picture($outdir);

exit;

# evalue classification
sub evalue_func{
        my ($e_v) = @_;
        if($e_v == 0){
                $ehash{'0'} += 1;
        }
        elsif($e_v =~ /(\d)+e\-(\d+)/){
                if($2 > 100){
                        $ehash{'0~1e-100'} += 1;
                }
                elsif($2 > 60){
                        $ehash{'1e-100~1e-60'} += 1;
                }
                elsif($2 > 45){
                        $ehash{'1e-60~1e-45'} += 1;
                }
                elsif($2 > 30){
                        $ehash{'1e-45~1e-30'} += 1;
                }
                elsif($2 > 15){
                        $ehash{'1e-30~1e-15'} += 1;
                }
                elsif($2 > 5){
                        $ehash{'1e-15~1e-5'} += 1;
                }
        }
        else{
                print "$e_v\n";
        }
}

# identity classification
sub similarity_func{
        my ($iden) = @_;
        if($iden > 0.95){
                $ihash{'95%~100%'} += 1;
        }
        elsif($iden > 0.8){
                $ihash{'80%~95%'} += 1;
        }
        elsif($iden > 0.6){
                $ihash{'60%~80%'} += 1;
        }
        elsif($iden > 0.4){
                $ihash{'40%~60%'} += 1;
        }
        else{
                $ihash{'0'} += 1;
        }
        if($iden < $min_iden){
                $min_iden = $iden;
        }
}


# obsolete function: specie classification
sub species_func{
        my ($annot) = @_;
        chomp $annot;
        if($annot =~ /.*\[\s*([^\[\]]+?)\s*\]/){
                $shash{$1} += 1;
        }
	else{
		$shash{'other'} += 1;
	}
}

# draw picture using plotrix
sub draw_picture{
        my ($odir) = @_;
        my $species_col1 = join("','", @species_array);
        $species_col1 = "'$species_col1'";
        my $species_col2;
        for (@species_array) { $species_col2 .= "$shash{$_}," }
        chop $species_col2; 
        open R,">$odir/draw_picture.R";
        print R<<RTXT;
library(ggplot2)

species = c($species_col1)
gene.numbers = c($species_col2)
percentage = paste(round(gene.numbers/sum(gene.numbers)*100, 2), '%', sep='')
species = data.frame(species, gene.numbers, percentage)
species[[1]] = factor(species[[1]], levels=species[[1]])

pdf("$odir/species_distribution.pdf", width=10, height=8)
ggplot(species, aes(x=factor(1), y=gene.numbers, fill=species))+
  geom_bar(stat='identity',width=1,,color='black')+
  geom_text(aes(x=1.75, y=cumsum(gene.numbers)-gene.numbers/2, label=percentage))+
  coord_polar(theta='y')+
  ggtitle('Species Distribution')+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        plot.title = element_text(face="bold",size=19))+
  scale_fill_manual(name='Species',values = c('royalblue','palegreen3','indianred2','deepskyblue','darkorange'))
dev.off()
RTXT
close R;
	system("Rscript $odir/draw_picture.R");
	system("convert -density 300 -resize 30% $odir/species_distribution.pdf $odir/species_distribution.png");
}
