#!/usr/bin/perl

=head1 Name

	auto_Secretory_proteinpl  -- the pipeline of Secretory_protein function annotation.

=head1 Description

	The path of the database and softwares invoked will be put in the "config.txt".

=head1 Version

	Author: lanlei (lanlei@aimigene.com).

=head1 Usage

	perl auto_Secretory_protein.pl [options] *.pep

	-workdir			set the workdir.
	-signalp			run Signalp.
	-tmhmm				run Tmhmm
	-type				Organism type> (euk, gram+, gram-). Default: 'euk'.
	-help				output help information to screen.

=head1 Exmple

	perl auto_Secretory_protein.pl -workdir ./ -signalp -tmhmm  pep.fasta
=cut

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin $Script);
use lib "$Bin/../../lib";
use GACP qw(parse_config);

my ($workdir,$type,,$signalp,$tmhmm,$help,$fasta);
GetOptions(
    "workdir:s"=>\$workdir,
	"type:s"=>\$type,
	"signalp"=>\$signalp,
	"tmhmm"=>\$tmhmm,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

$fasta = shift;
$workdir ||= "./";
$type ||= "gram+";

$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

my $time = `date`; chomp($time);
print "Start Secretory_protein function annotation at ... $time! \n";

my $config_file = "$Bin/../../config.txt";
my $Signalp_db = parse_config($config_file,"Signalp_db");
my $Tmhmm_soft = parse_config($config_file,"Tmhmm_soft");
my $tmp_dir = "$workdir/tmp";

if (defined $signalp){
	&mkdir_chdir("$workdir/Signalp");
	open SH, ">$workdir/Signalp/$filename.signalp.sh" || die $!;
	print SH "#!/bin/bash\n";
	print SH "export PATH=$Signalp_db:\$PATH \n";
	print SH "mkdir -p $tmp_dir \n";
	print SH "signalp -batch 30000 -org $type -fasta $fasta -tmp $tmp_dir -gff3 -mature -format long \n";

	system("sh $workdir/Signalp/$filename.signalp.sh 1>$workdir/Signalp/$filename.signalp.sh.e 2>$workdir/Signalp/$filename.signalp.sh.o ") == 0 || die $!;
}

if (defined $tmhmm){
	&mkdir_chdir("$workdir/Tmhmm");
	open SH, ">$filename.tmhmm.sh" || die $!;
	print SH "#!/bin/bash\n";
	print SH "$Tmhmm_soft $fasta > $workdir/Tmhmm/tmm.txt \n";

	system("sh $workdir/Tmhmm/$filename.tmhmm.sh 1>$workdir/Tmhmm/$filename.tmhmm.sh.e 2>$workdir/Tmhmm/$filename.tmhmm.sh.o ") == 0 || die $!;
}

print "Secretory_protein annotation DONE ... $time!\n";

sub mkdir_chdir{
	my $dir_name = shift;
	system("rm -rf $dir_name") == 0 || die $! if (-d $dir_name);
	system("mkdir $dir_name") == 0 || die $!;
	chdir $dir_name or die ("ERROR: Cannot change directory to $dir_name: $!\n");
}