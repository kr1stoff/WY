#!/usr/bin/perl

=head1 Name

	auto_Pfampl  -- the pipeline of Pfam function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_Pfam.pl [options] *.pep

	-workdir			set the workdir. 
	-help				output help information to screen.

=head1 Exmple

	perl auto_Pfam.pl -workdir ./ pep.fasta
=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

my ($workdir,$help,$fasta);
GetOptions(
    "workdir:s"=>\$workdir,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

$fasta = shift;
$workdir ||= "./";

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

#&mkdir_chdir("$workdir/Pfam");
my $time = `date`; chomp($time);
print "Start Pfam function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $lib = parse_config($config_file,"lib");
my $Pfam_db = parse_config($config_file,"Pfam_db");

open SH, ">$workdir/$filename.Pfam.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "export PERL5LIB=$lib/:\${PERL5LIB}\n";
print SH "perl $func_bin/pfam_scan.pl -fasta $fasta -dir $Pfam_db -outfile $filename.pfam.out \n";
close SH;

system("sh $filename.Pfam.sh") == 0 || die $!;
print "Pfam annotation DONE ... $time!\n";
