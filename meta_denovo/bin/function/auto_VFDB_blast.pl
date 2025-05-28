#!/usr/bin/perl

=head1 Name

	auto_VFDB_blast.pl  -- the pipeline of VFDB function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_VFDB_blast.pl [options] *.pep

	-Evalue				set the blast evalue, default=1e-5.
	-workdir			set the workdir.
	-tools				set the mapping tools(diamond|blastp),default="blastp".
	-help				output help information to screen.

=head1 Exmple

	perl auto_VFDB_blast.pl -tools blastp -$workdir ./ -Evalue 1e-5 pep.fasta
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

#&mkdir_chdir("$workdir/VFDB");
my $time = `date`; chomp($time);
print "Start VFDB function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $VFDB_db = parse_config($config_file,"VFDB");
my $VFDB_anno = parse_config($config_file,"VFDB_anno");

# my $outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname qlen slen'";

open SH, ">$workdir/$filename.VFDB.sh" || die $!;
print SH "#!/bin/bash\n";

if ($tools eq "diamond"){
	print SH "diamond blastp -q $fasta -d $VFDB_db --quiet -p 10 -e $Evalue -k 50 --sensitive -o $filename.blast.VFDB.out\n";
}
if ($tools eq "blastp"){
	print SH "blastp -db $VFDB_db -query $fasta -evalue $Evalue -num_threads 10 -outfmt 6 -out $filename.blast.VFDB.out\n";
}

print SH "perl $func_bin/0.choose_blast_m8.pl -i $filename.blast.VFDB.out -o $filename.blast.VFDB.out.filter -b Y -p 40 -s $VFDB_db -q $fasta \n";
print SH "perl $func_bin/1.anno.pl $filename.blast.VFDB.out.filter $filename.blast.VFDB.out.filter.xls $VFDB_anno \n";
close SH;

system("sh $filename.VFDB.sh") == 0 || die $!;
print "VFDB annotation DONE ... $time!\n";
