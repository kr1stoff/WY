#!/usr/bin/perl

=head1 Name

	auto_CARD.pl  -- the pipeline of CARD function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei.

=head1 Usage

	perl auto_CARD.pl [options] *.pep

	-workdir			set the workdir.
	-help				output help information to screen.

=head1 Exmple

	perl auto_CARD.pl -workdir ./ pep.fasta
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

#&mkdir_chdir("$workdir/CARD");
my $time = `date`; chomp($time);
print "Start CARD function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $CARD_db = parse_config($config_file,"CARD_db");
my $activate = parse_config($config_file,"activate");


open SH, ">$workdir/$filename.CARD.sh" || die $!;
print SH "#!/bin/bash \n";
print SH "source $activate rgi \n";
print SH "rgi load --card_json $CARD_db/card.json --local\n";
print SH "rgi main --input_sequence $fasta --output_file $workdir/out --input_type protein --local --clean \n";
close SH;

system("bash $filename.CARD.sh 1>$filename.CARD.sh.e 2>$filename.CARD.sh.o");
print "CARD annotation DONE ... $time!\n";
