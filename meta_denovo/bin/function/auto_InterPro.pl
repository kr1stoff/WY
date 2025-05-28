#!/usr/bin/perl

=head1 Name

	auto_InterPropl  -- the pipeline of InterPro function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_InterPro.pl [options] *.pep

	-workdir			set the workdir.
	-appl				set applications, this option can be used multiple times
	-help				output help information to screen.

=head1 Exmple

	perl auto_InterPro.pl -workdir ./ pep.fasta
=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

my ($workdir,$appl,$prefix,$help,$fasta);
GetOptions(
    "workdir:s"=>\$workdir,
	"appl:s"=>\$appl,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

$fasta = shift;
$workdir ||= "./";
$appl ||= "--appl ProDom --appl PRINTS --appl Pfam --appl SMART --appl PANTHER --appl ProSiteProfiles --appl ProSitePatterns ";

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

#&mkdir_chdir("$workdir/InterPro");
my $time = `date`; chomp($time);
print "Start InterPro function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $InterPro_script = parse_config($config_file,"InterPro_script");

open SH, ">$workdir/$filename.InterPro.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "source /share/database/.config/.bash_profile\n";
print SH "$InterPro_script -i $fasta -cpu 10 -f tsv -goterms $appl -o $workdir/$filename.iprscan\n";
print SH "perl $func_bin/iprscan_parser_xls.pl $workdir/$filename.iprscan $workdir/$filename.iprscan.xls\n";
print SH "perl $func_bin/iprscan_parser51-55.pl $workdir/$filename.iprscan -outdir $workdir\n";
close SH;

system("sh $filename.InterPro.sh") == 0 || die $!;
print "InterPro annotation DONE ... $time!\n";
