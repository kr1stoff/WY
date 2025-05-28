#!/usr/bin/perl

=head1 Name

	auto_PHI_blast.pl  -- the pipeline of PHI function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_PHI_blast.pl [options] *.pep

	-Evalue				set the blast evalue, default=1e-5.
	-workdir			set the workdir. 
	-tools				set the mapping tools(diamond|blastp),default="blastp".
	-help				output help information to screen.

=head1 Exmple

	perl auto_PHI_blast.pl -workdir ./ -tools blastp -Evalue 1e-5 pep.fasta
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

#&mkdir_chdir("$workdir/PHI");
my $time = `date`; chomp($time);
print "Start PHI function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $PHI_db = parse_config($config_file,"PHI");
my $PHI_anno = parse_config($config_file,"PHI_anno");

# my $outfmt = "'6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname qlen slen'";

open SH, ">$workdir/$filename.PHI.sh" || die $!;
print SH "#!/bin/bash\n";

if ($tools eq "diamond"){
	print SH "diamond blastp -q $fasta -d $PHI_db --quiet -p 10 -e $Evalue -k 50 --sensitive -o $filename.blast.PHI.out\n";
}
if ($tools eq "blastp"){
	print SH "blastp -db $PHI_db -query $fasta -evalue $Evalue -num_threads 10 -outfmt 6 -out $filename.blast.PHI.out\n";
}

print SH "perl $func_bin/0.choose_blast_m8.pl -i $filename.blast.PHI.out -o $filename.blast.PHI.out.filter -b Y -p 40 -s $PHI_db -q $fasta \n";
print SH "perl $func_bin/1.anno.pl $filename.blast.PHI.out.filter $filename.blast.out.filter.xls $PHI_anno \n";
close SH;

system("sh $filename.PHI.sh") == 0 || die $!;
print "PHI annotation DONE ... $time!\n";
