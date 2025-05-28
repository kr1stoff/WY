#!/usr/bin/perl
=head1 Name
	auto_fun_ann.pl  -- the pipeline of function annotation.

=head1 Description
	This program invoke GO,NR2GO,KEGG,COG,eggNOG,CAZY,,VFDB,CARD,PHIbase,TCDB..
	The path of the softwares invoked will be put in the "config.txt".

=head1 Usage
	perl auto_fun_ann.pl [options] *.pep
	note: the pep_file must be input by absolute path

    -NR             run NR.
    -GO             run GO.
      -NR2GO        run GO based on NR.
      -eggNOG2GO    run GO based on eggNOG.
    -KEGG           run KEGG.
    -COG            run COG.
    -eggNOG         run eggNOG.
    -CAZY           run CAZY.
    -VFDB           run VFDB.
    -CARD           run CARD.
    -PHIbase        run PHIbase.
    -TCDB           run TCDB.
    -species        set the species type, bacteria, animal or plant; if not set, will compare to bacteria type.
    -Evalue         set the blast evalue, default=1e-5.
    -help           output help information to screen.

=head1 Exmple
nohup perl auto_fun_ann.pl -Interpro -KEGG -Swissprot -Trembl -KOG -NT -NR -GO -species bacteria pep.fasta 

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/../lib";
use GACP qw(parse_config);
use Cwd qw(getcwd abs_path);

my ($NR,$GO,$NR2GO,$eggNOG2GO,$KEGG,$COG,$eggNOG,$CAZY,$VFDB,$CARD,$PHIbase,$TCDB);
my ($species, $evalue,$workdir);
my ($help);

GetOptions(
	"NR" =>\$NR,
	"GO" =>\$GO,
		"NR2GO" =>\$NR2GO,
		"eggNOG2GO" =>\$eggNOG2GO,
	"KEGG"=>\$KEGG,
	"COG" =>\$COG,
	"eggNOG" =>\$eggNOG,
	"CAZY" =>\$CAZY,
	"VFDB" =>\$VFDB,
	"CARD" =>\$CARD,
	"PHIbase" =>\$PHIbase,
	"TCDB" =>\$TCDB,
	"species:s"=>\$species,
	"Evalue"=>\$evalue,
	"workdir"=>\$workdir,
	"Help"=>\$help,
);

die `pod2text $0` if(@ARGV==0 || $help);

$evalue ||= 1e-5;
$species ||="Bacteria";
$workdir ||="./";

my $dir = `pwd`; chomp ($dir);
my $fasta = shift;
$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $fasta_name = basename($fasta);
my $scomnames = $fasta_name;
my $finish_string = "All jobs finish";
my $config_file = "$Bin/../config.txt";

my $func_bin = parse_config($config_file,"func_bin");
my $GO_path = parse_config($config_file,"GO");

my $NR_blast_tools ||= "diamond";

open(OUT , ">" . "$workdir/STEP01_fun_ann_work.sh") or die $!;

# if (defined $TCDB)
# {
# 	print OUT "########## TCDB ##########\n";
# 	print OUT "mkdir TCDB\ncd TCDB\n";
# 	print OUT "nohup perl $func_bin/auto_TCDB_blast.pl -Evalue 1e-5 -tools blastp $fasta && \\\n";
# 	print OUT "cd ..\n";
# }

# if (defined $VFDB)
# {
# 	print OUT "########## VFDB ##########\n";
# 	print OUT "mkdir VFDB\ncd VFDB\n";
# 	print OUT "nohup perl $func_bin/auto_VFDB_blast.pl -Evalue 1e-5 -tools blastp $fasta && \\\n";
# 	print OUT "cd ..\n";
# }

# if (defined $PHIbase)
# {
# 	print OUT "########## PHIbase ##########\n";
# 	print OUT "mkdir PHIbase\ncd PHIbase\n";
# 	print OUT "nohup perl $func_bin/auto_PHI_blast.pl -Evalue 1e-5 -tools blastp $fasta && \\\n";
# 	print OUT "cd ..\n";
# }

if (defined $CARD)
{
	print OUT "########## CARD ##########\n";
	print OUT "mkdir -p $workdir/CARD\ncd $workdir/CARD\n";
	print OUT "nohup perl $func_bin/auto_CARD.pl $fasta && \\\n";
	print OUT "cd ..\n";
}

# if (defined $CAZY)
# {
# 	print OUT "########## CAZY ##########\n";
# 	print OUT "mkdir -p $workdir/CAZY\ncd $workdir/CAZY\n";
# 	print OUT "nohup perl $func_bin/auto_CAZY.pl -workdir ./ $fasta && \\\n";
# 	print OUT "cd ..\n";
# }

if (defined $KEGG || defined $COG || $eggNOG2GO)
{
	print OUT "########## eggNOG ##########\n";
	print OUT "nohup perl $func_bin/auto_eggNOG.pl -eggNOG $fasta -workdir $workdir && \\\n";
}

if (defined $eggNOG && !defined $KEGG && !defined $COG && !defined $eggNOG2GO)
{
	print OUT "########## eggNOG ##########\n";
	print OUT "nohup perl $func_bin/auto_eggNOG.pl -eggNOG $fasta -workdir $workdir && \\\n";
}

if (defined $KEGG)
{
	print OUT "########## KEGG ##########\n";
	print OUT "nohup perl $func_bin/auto_eggNOG.pl -KEGG $fasta -workdir $workdir && \\\n";
}

if (defined $COG)
{
	print OUT "########## COG ##########\n";
	print OUT "nohup perl $func_bin/auto_eggNOG.pl -COG $fasta -workdir $workdir && \\\n";
}

if (defined $GO)
{
	if (defined $NR2GO)
	{
		print OUT "########## NR && GO ##########\n";
		print OUT "nohup perl $func_bin/auto_nr_blast3GO.pl $fasta $species $evalue ./ $NR_blast_tools GO && \\\n";
	}
	if (defined $eggNOG2GO)
	{
		print OUT "########## GO ##########\n";
		print OUT "nohup perl $func_bin/auto_eggNOG.pl -GO $fasta -workdir $workdir && \\\n";
	}
}

if (defined $NR && ! defined $GO)
{
	print OUT "########## NR ##########\n";
	print OUT "nohup perl $func_bin/auto_nr_blast3GO.pl $fasta $species 1e-5 ./ $NR_blast_tools NO ./ && \\\n";
}

print OUT "echo $finish_string > STEP01_fun_ann_work.sh.sign\n";

close(OUT);

open(OUT , ">" . "$workdir/STEP02_fun_ann_stat.sh") or die $!;
	my $stat_parameter = "";
	print OUT "cd $dir\n";
	print OUT "mkdir -p $workdir/fun_ann_stat\n";
	print OUT "cd $workdir/fun_ann_stat\n";
	print OUT "ln -s $fasta $fasta_name && \\\n";
	print OUT "grep \'^>\' $fasta_name | sed \'s/>//\' | awk \'{print \$1}\' > all_gene.id && \\\n";
	if (defined $NR)
	{
		print OUT "ln -s $dir/NR/$scomnames.blast.nr.xls && \\\n";
		$stat_parameter .= "-nr $scomnames.blast.nr.xls ";
	}
	if (defined $GO)
	{
		if (defined $NR2GO)
		{
			print OUT "ln -s  $dir/GO_NR/$scomnames.Gene2GO.xls && \\\n";
			$stat_parameter .= "-go $scomnames.Gene2GO.xls -obo $GO_path/go-basic.obo ";
		}
		if (defined $eggNOG2GO)
		{
			print OUT "ln -s  $dir/GO_eggNOG/$scomnames.Gene2GO.xls && \\\n";
			$stat_parameter .= "-go $scomnames.Gene2GO.xls -obo $GO_path/db/go-basic.obo ";
		}
	}
	if (defined $eggNOG)
	{
			print OUT "ln -s  $dir/eggNOG/out.emapper.annotations && \\\n";
			$stat_parameter .= "-eggNOG out.emapper.annotations ";
	}
	if (defined $KEGG)
	{
		print OUT "ln -s $dir/KEGG/$scomnames.kegg.xls && \\\n";
		$stat_parameter .= "-kegg $scomnames.kegg.xls ";
	}


	print OUT "perl $func_bin/all_function_stat.pl -list all_gene.id $stat_parameter -outxls ./annotation.xls --outstat ./annotation_stat.xls &&\\\n";
	print OUT "cd ..\n";
	print OUT "echo $finish_string > STEP02_fun_ann_stat.sh.sign\n";

close(OUT);

# open(OUT , ">" . "$workdir/STEP03_delete_tmp_files.sh") or die $!;
# #	print OUT "cd function_annotation\n";
# print OUT "echo $finish_string > STEP03_delete_tmp_files.sh.sign\n";
# close(OUT);
