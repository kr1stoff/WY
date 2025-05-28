#!/usr/bin/perl

=head1 Name

	auto_ARDB_blast.pl  -- the pipeline of ARDB function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_ARDB_blast.pl [options] *.pep

	-workdir			set the workdir. 
	-prefix				set the prefix of the output prefix.
	-help				output help information to screen.

=head1 Exmple

	perl auto_ARDB_blast.pl -workdir ./ -prefix test pep.fasta
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
	"prefix:s"=>\$prefix,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

$fasta = shift;
$workdir ||= "./";

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

#&mkdir_chdir("$workdir/ARDB");
my $time = `date`; chomp($time);
print "Start ARDB function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $ARDB_db = parse_config($config_file,"ARDB_db");

open INFO, ">$workdir/genomeList.tab";
print INFO ">$prefix\n$fasta\n";
close INFO;

open SH, ">$workdir/$filename.ARDB.sh" || die $!;
print SH "#!/bin/bash\n";
print SH "export PERL5LIB=$ARDB_db/:\${PERL5LIB} \n";
print SH "perl $func_bin/ardbAnno_modified.pl genomeList.tab \n";
close SH;

system("sh $filename.ARDB.sh") == 0 || die $!;
print "ARDB annotation DONE ... $time!\n";
