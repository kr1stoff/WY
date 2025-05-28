#!/usr/bin/perl

=head1 Name

	auto_TCDB_blast.pl  -- the pipeline of TCDB function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_TCDB_blast.pl [options] *.pep

	-Evalue				set the blast evalue, default=1e-5.
	-workdir			set the workdir. 
	-tools				set the mapping tools(diamond|blastp),default="blastp".
	-help				output help information to screen.

=head1 Exmple

	perl auto_TCDB_blast.pl -tools blastp -workdir ./ -Evalue 1e-5 pep.fasta
=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

my ($Evalue,$workdir,$tools,$help,$fasta);
GetOptions(
    "Evalue:s"=>\$Evalue,
    "workdir:s"=>\$workdir,
	"tools:s"=>\$tools,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

$fasta = shift;
$Evalue ||= "1e-5";
$workdir ||= "./";
$tools ||= "blastp";

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

#&mkdir_chdir("$workdir/TCDB");
my $time = `date`; chomp($time);
print "Start TCDB function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $TCDB_db = parse_config($config_file,"TCDB");
my $TCDB_anno = parse_config($config_file,"TCDB_anno");

# my $outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname qlen slen'";

open SH, ">$workdir/$filename.TCDB.sh" || die $!;
print SH "#!/bin/bash\n";

if ($tools eq "diamond"){
	print SH "diamond blastp -q $fasta -d $TCDB_db --quiet -p 10 -e $Evalue -k 50 --sensitive -o $filename.blast.TCDB.out\n";
}
if ($tools eq "blastp"){
	print SH "blastp -db $TCDB_db -query $fasta -evalue $Evalue -num_threads 10 -outfmt 6 -out $filename.blast.TCDB.out\n";
}

print SH "perl $func_bin/0.choose_blast_m8.pl -i $filename.blast.TCDB.out -o $filename.blast.TCDB.out.filter -b Y -p 40 -s $TCDB_db -q $fasta \n";
# add TCDB class information
print SH "perl $func_bin/1.anno_TCDB.pl $filename.blast.TCDB.out.filter $filename.blast.TCDB.out.filter.xls $TCDB_anno \n";
close SH;

system("sh $filename.TCDB.sh") == 0 || die $!;
print "TCDB annotation DONE ... $time!\n";
