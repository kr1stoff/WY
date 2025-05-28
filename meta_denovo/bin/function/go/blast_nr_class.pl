#!/usr/bin/perl -w

=head1 Name

	blast_nr_class.pl  --  statistics of blastx nr result

=head1 Description
  This program is designed for doing statistics of blastx nr result.

=head1 Version

  Vertion:  3.0
  Date:	    2013-07-22
  Author:   Xu Tong
=head1 Usage

  perl blast_nr_class.pl <options>

=head1 sample 
	perl blast_nr_class.pl -nr xxxx.fa.blast.Nr.xls -nohead -out ./
=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename;
use DBI;
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
#	print STDERR "$total_gene\n";
        my @nr_arr = split(/\t/, $_,13);
        ($gene_name, $nr_id, $identity, $evalue, $annotation) = ($nr_arr[0], $nr_arr[1], $nr_arr[2], $nr_arr[10], $nr_arr[12]);
#        evalue_func($evalue);
#		my $identity_ptg = $identity/100;
#        similarity_func($identity_ptg);
#	species_func_db($nr_id);	
	species_func($annotation);
}
close NR;

=modify
# output the evalue info
open EOUT,">$outdir/evalue_statistic.xls";
print EOUT "evalue\tgene numbers\tpercentage\n";
my @evalue_array = ('0', '0~1e-100', '1e-100~1e-60', '1e-60~1e-45', '1e-45~1e-30', '1e-30~1e-15', '1e-15~1e-5');
foreach my $eKey(@evalue_array){
	my $per = 0;
	if($ehash{$eKey}){
	        $per = 100 * $ehash{$eKey} / $total_gene;
	}else{
		$ehash{$eKey} = 0;
	}
        printf EOUT "%s\t%d\t%.2f%s\n", $eKey, $ehash{$eKey}, $per, "%";
}
close EOUT;

# output the identity info
open IOUT,">$outdir/similarity_statistic.xls";
print IOUT "similarity\tgene numbers\tpercentage\n";
foreach my $iKey(sort keys %ihash){
        my $iPer = 100 * $ihash{$iKey} / $total_gene;
        if($iKey eq '0'){
                $min_iden = int($min_iden * 100);
                printf IOUT "%s%s\t%d\t%.2f%s\n", $min_iden, "%~40%", $ihash{$iKey}, $iPer, "%";
        }
        else{
                printf IOUT "%s\t%d\t%.2f%s\n", $iKey, $ihash{$iKey}, $iPer, "%";
        }
}
close IOUT;
=cut
# output the species info
our @species_array;
my $species_count = 0;
#my $sOtherPer = 100;
open SOUT,">$outdir/species_statistic.xls";
print SOUT "species\tgene numbers\tpercentage\n";
foreach my $sKey(sort {$shash{$b} <=> $shash{$a}} keys %shash){
        my $sPer = sprintf("%.2f", 100 * $shash{$sKey} / $total_gene);
	if($sKey ne "other"){
#		$species_count += 1;
		printf SOUT "%s\t%d\t%.2f%s\n", $sKey, $shash{$sKey}, $sPer, "%";
		if($species_count < 4 && $sPer >= 1){
			$species_count += 1;
#        	        $sOtherPer -= $sPer;
#                	printf SOUT "%s\t%d\t%.2f%s\n", $sKey, $shash{$sKey}, $sPer, "%";
	                push @species_array, $sKey;
        	}
	        else{
        	        $shash{'other'} += $shash{$sKey};
		}
        }
}
#print STDERR $shash{'other'};
push @species_array, 'other';
#printf SOUT "%s\t%d\t%.2f%s\n", "other", $shash{'other'}, $sOtherPer, "%";

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

# species classification using DBI
=del
sub species_func_db{
	my ($gi) = @_;
	$gi =~ s/^gi\|(\d+)\|.*$/$1/g;
	my $table = "gi_taxid_name_".sprintf("%03d", int($gi / 2000000));
#	print STDERR "select * from $table where gi=$gi\n";
	my $sql = "select * from $table where gi=$gi";	
#	my $sql = "select * from $table left join names on $table.tax_id=names.tax_id where $table.gi=$gi";
	my $sth = $dbh->prepare($sql);
	$sth->execute() or die $!;
	my ($taxid, $taxname);
	$sth->bind_columns(undef, \$gi, \$taxid, \$taxname);
	$sth->fetch();
	$sth->finish();
	if($taxname eq "Ricimus communis"){
		$taxname = "Ricinus communis";
	}
	$shash{$taxname} += 1;
}
=cut
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

#evalue = read.table('$odir/evalue_statistic.xls', sep='\t', header=T)
#similarity = read.table('$odir/similarity_statistic.xls', sep='\t', header=T)
species = c($species_col1)
gene.numbers = c($species_col2)
percentage = paste(round(gene.numbers/sum(gene.numbers)*100, 2), '%', sep='')
species = data.frame(species, gene.numbers, percentage)
species[[1]] = factor(species[[1]], levels=species[[1]])

#pdf("$odir/evalue_distribution.pdf", width=10, height=10)
#ggplot(evalue, aes(x=factor(1), y=gene.numbers, fill=evalue))+
#  geom_bar(stat='identity',width=1,color='black')+
#  geom_text(aes(x=1.75, y=cumsum(gene.numbers)-gene.numbers/2, label=percentage))+
#  coord_polar(theta='y')+
#  ggtitle('E-value Distribution')+
#  theme(panel.grid=element_blank(),
#        panel.background=element_blank(),
#        axis.title=element_blank(),
#        axis.text=element_blank(),
#        axis.ticks=element_blank(),
#        legend.title=element_text(size=16),
#        legend.text=element_text(size=16),
#        plot.title = element_text(face="bold"))+
#  scale_fill_discrete(name='E-value')
#dev.off()

#pdf("$odir/similarity_distribution.pdf", width=10, height=10)
#ggplot(similarity, aes(x=factor(1), y=gene.numbers, fill=similarity))+
#  geom_bar(stat='identity',width=1,color='black')+
#  geom_text(aes(x=1.75, y=cumsum(gene.numbers)-gene.numbers/2, label=percentage))+
#  coord_polar(theta='y')+
#  ggtitle('Similarity Distribution')+
#  theme(panel.grid=element_blank(),
#        panel.background=element_blank(),
#        axis.title=element_blank(),
#        axis.text=element_blank(),
#        axis.ticks=element_blank(),
#        legend.title=element_text(size=16),
#        legend.text=element_text(size=16),
#        plot.title = element_text(face="bold"))+
#  scale_fill_discrete(name='Similarity')
#dev.off()

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
#	system("$Bin/../software/convert -density 300 -resize 40% $odir/evalue_distribution.pdf $odir/evalue_distribution.png");
#	system("$Bin/../software/convert -density 300 -resize 40% $odir/similarity_distribution.pdf $odir/similarity_distribution.png");
	system("convert -density 300 -resize 30% $odir/species_distribution.pdf $odir/species_distribution.png");
}
