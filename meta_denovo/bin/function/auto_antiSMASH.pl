#!/usr/bin/perl

=head1 Name

	auto_antiSMASH.pl  -- the pipeline of antiSMASH function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_antiSMASH.pl [options] *.pep

	-workdir			set the workdir.
	-help				output help information to screen.

=head1 Exmple

	perl auto_antiSMASH.pl -workdir ./ pep.fasta
=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

my ($workdir,$prefix,$help,$gbk);
GetOptions(
    "workdir:s"=>\$workdir,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

$gbk = shift;
$workdir ||= "./";

$gbk = abs_path($gbk);
$workdir = abs_path($workdir);
my $filename = basename($gbk);

#&mkdir_chdir("$workdir/antiSMASH");
my $time = `date`; chomp($time);
print "Start antiSMASH function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";

open SH, ">$workdir/$filename.antiSMASH.sh" || die $!;
print SH "#!/bin/bash \n";
print SH "source activate antismash \n";
print SH "antismash --clusterhmmer --cb-subcluster --cb-knownclusters --cb-general --asf --pfam2go --smcog-trees --cpus 10 $gbk \n";
print SH "conda deactivate \n";
close SH;

system("sh $filename.antiSMASH.sh 1>$filename.antiSMASH.sh.e 2>$filename.antiSMASH.sh.o") == 0 || die $!;
print "antiSMASH annotation DONE ... $time!\n";
