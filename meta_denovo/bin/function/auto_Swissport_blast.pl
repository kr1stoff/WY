#!/usr/bin/perl

=head1 Name

	auto_Swissport_blast.pl  -- the pipeline of Swissport function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_Swissport_blast.pl [options] *.pep

	-dbClass			set the species type, Bacteria or Viruses; if not set, will compare to all type.
	-Evalue				set the blast evalue, default=1e-5.
	-workdir			set the workdir.
	-tools				set the mapping tools(diamond|blastp),default="blastp".
	-help				output help information to screen.

=head1 Exmple

	perl auto_Swissport_blast.pl -dbClass Bacteria -Evalue 1e-5 -tool blastp -workdir ./ spep.fasta
=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

my ($dbClass,$Evalue,$workdir,$tools,$help,$fasta);
GetOptions(
    "dbClass:s" =>\$dbClass,
    "Evalue:s"=>\$Evalue,
    "workdir:s"=>\$workdir,
	"tools:s"=>\$tools,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

$fasta = shift;
$dbClass ||= "all";
$Evalue ||= "1e-5";
$workdir ||= "./";
$tools ||= "blastp";

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

#&mkdir_chdir("$workdir/Swissport");
my $time = `date`; chomp($time);
print "Start Swissport function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $Swissport = parse_config($config_file,"Swissport");
my $Swissport_anno = parse_config($config_file,"Swissport_anno");

my $Swissport_db;
if ($dbClass eq "Bacteria"){ $Swissport_db = "$Swissport/uniprot_sprot.Bacteria.fasta.simple"; }
if ($dbClass eq "Viruses"){ $Swissport_db = "$Swissport/uniprot_sprot.Viruses.fasta.simple"; }
if ($dbClass eq "all"){ $Swissport_db = "$Swissport/uniprot_sprot.fasta.simple"; }

# my $outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname qlen slen'";

open SH, ">$workdir/$filename.Swissport.sh" || die $!;
print SH "#!/bin/bash\n";

if ($tools eq "diamond"){
	print SH "diamond blastp -q $fasta -d $Swissport_db --quiet -p 10 -e $Evalue -k 50 --sensitive -o $filename.blast.swissprot.out\n";
}
if ($tools eq "blastp"){
	print SH "blastp -db $Swissport_db -query $fasta -evalue $Evalue -num_threads 10 -outfmt 6 -out $filename.blast.swissprot.out\n";
}

print SH "perl $func_bin/Swissport_get_anno.pl -tophit 1 -topmatch 1 -id $Swissport_anno -input $filename.blast.swissprot.out -out $filename.blast.swissprot.xls \n";
close SH;

system("sh $filename.Swissport.sh") == 0 || die $!;
print "Swissport annotation DONE ... $time!\n";
