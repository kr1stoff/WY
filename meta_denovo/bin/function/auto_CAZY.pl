#!/usr/bin/perl

=head1 Name

	auto_CAZYpl  -- the pipeline of CAZY function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_CAZY.pl [options] *.pep

	-workdir			set the workdir.
	-help				output help information to screen.

=head1 Exmple

	perl auto_CAZY.pl -workdir ./ pep.fasta
=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

my ($workdir,$prefix,$help,$fasta);
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

#&mkdir_chdir("$workdir/CAZY");
my $time = `date`; chomp($time);
print "Start CAZY function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $CAZY_db = parse_config($config_file,"CAZY_db");

open SH, ">$workdir/$filename.CAZY.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "source activate run_dbcan\n";
print SH "run_dbcan.py --out_dir $workdir/out_$filename --db_dir $CAZY_db $fasta protein \n";
close SH;

system("sh $filename.CAZY.sh 1>$filename.CAZY.sh.e 2>$filename.CAZY.sh.o") == 0 || die $!;
print "CAZY annotation DONE ... $time!\n";
