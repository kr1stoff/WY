#!/usr/bin/perl

use warnings;
use strict;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

die "perl $0 <gene.fasta> <Bacteria|Fungi|Viruses> <evalue> <workdir> <diamond|blastp> <GO|NO>" if (@ARGV != 6);
my $fasta = shift;
my $dbClass = shift;
my $Evalue = shift;
my $workdir = shift;
my $tools = shift;
my $GO = shift;

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

&mkdir_chdir("$workdir/NR");
my $time = `date`; chomp($time);
print "Start NR function annotation at ... $time!\nMake directory  $workdir/NR && cd  $workdir/NR ...\n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $go_bin = "$Bin/go";
my $BLASTDB=parse_config($config_file,"BLASTDB");
my $NR_path = parse_config($config_file,"NR");
my $GO_path = parse_config($config_file,"GO");

my $nr_db;
if ($dbClass eq "Bacteria"){ $nr_db = "$NR_path/Bacteria.fa"; }
if ($dbClass eq "Fungi"){ $nr_db = "$NR_path/Fungi.fa"; }
if ($dbClass eq "Viruses"){$nr_db = "$NR_path/Viruses.fa"; }

# my $outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname qlen slen'";

open SH_blast, ">$workdir/NR/$filename.blast.sh" || die $!;
if ($tools eq "diamond"){
	print SH_blast "diamond blastp -q $fasta -d $nr_db -p 20 -e $Evalue -k 50 --sensitive --quiet -o $workdir/NR/$filename.blast.nr \n";
}
if ($tools eq "blastp"){
	print SH_blast "export BLASTDB=$BLASTDB\n";
	print SH_blast "blastp -db $nr_db -query $fasta -evalue $Evalue -num_threads 20 -outfmt 6 -out $workdir/NR/$filename.blast.nr \n";
}
close SH_blast;

system("sh $workdir/NR/$filename.blast.sh") == 0 || die $!;

open SH_NR, ">$workdir/NR/blast.nr.process.sh" || die $!;
print SH_NR "#!/bin/bash\n";
print SH_NR "perl $func_bin/get_annot.pl -tophit 1 -topmatch 1 -id $nr_db.anno -input $workdir/NR/$filename.blast.nr -out $workdir/NR/$filename.blast.nr.xls && \\\n";
print SH_NR "perl $func_bin/blast_nr_class.pl -nr $workdir/NR/$filename.blast.nr.xls -outdir $workdir/NR \n";
close SH_NR;

system("sh $workdir/NR/blast.nr.process.sh") == 0 || die $!;

$time = `date`; chomp($time);
if ($GO =~ /GO/i){
	&mkdir_chdir("$workdir/GO_NR");
	print "Start GO function annotation at $time!\nMake directory  $workdir/GO_NR && cd  $workdir/GO_NR ..\n";

	open SH_GO, ">$workdir/GO_NR/blast3go.process.sh" || die $!;
	print SH_GO "#!/bin/bash\n";
	print SH_GO "ln -s $workdir/NR/$filename.blast.nr $workdir/GO_NR/$filename.blast.nr.annot && \\\n";
	print SH_GO "ln -s $workdir/NR/$filename.blast.nr.xls $workdir/GO_NR/$filename.blast.nr.annot.xls && \\\n";
	print SH_GO "python $go_bin/BLAST3GO.py --NR_tab $workdir/GO_NR/$filename.blast.nr.annot.xls --Accession_go $GO_path/db/accession2go --output $workdir/GO_NR/$filename.blast3go.annot && \\\n";
	print SH_GO "if [ -d $workdir/GO_NR/data ];then rm -rf $workdir/GO_NR/data;fi && \\\n";
	print SH_GO "mkdir -p $workdir/GO_NR/data && \\\n";
	print SH_GO "perl $go_bin/annot2goa.pl $GO_path/db/go-basic.obo $workdir/GO_NR/$filename.blast3go.annot $workdir/GO_NR/data/species && \\\n";
	print SH_GO "cat $workdir/GO_NR/data/species.[CFP] | cut -f 2 | sort -u > $workdir/GO_NR/gene.list && \\\n";
	print SH_GO "perl $go_bin/drawGO.pl -list $workdir/GO_NR/gene.list -goprefix  $workdir/GO_NR/data/species -goclass $GO_path/db/go.class -outprefix $workdir/GO_NR/$filename\n";
	close SH_GO;
	system("sh $workdir/GO_NR/blast3go.process.sh") == 0 || die $!;
	$time = `date`; chomp($time);
	print "NR && Go annotation DONE ... $time!\n";
	
}else{
	print "NR annotation DONE ... $time!\n";
}

#####################################sub#################################
sub mkdir_chdir{
	my $dir_name = shift;
	system("rm -rf $dir_name") == 0 || die $! if (-d $dir_name);
	system("mkdir $dir_name") == 0 || die $!;
	chdir $dir_name or die ("ERROR: Cannot change directory to $dir_name: $!\n");
}
