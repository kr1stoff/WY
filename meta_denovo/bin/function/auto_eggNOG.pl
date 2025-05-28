#!/usr/bin/perl
=head1 Name
	auto_eggNOG.pl  -- the pipeline of function annotation.

=head1 Description
	auto run the eggNOG annotation pipeline

=head1 Usage
	perl auto_eggNOG.pl [options] *.pep
	note: the pep_file must be input by absolute path

	-workdir		the work directory
	-eggNOG			run eggNOG
	-GO  			run GO
	-KEGG			run KEGG
	-COG 			run COG
	-help			output help information to screen

=head1 Exmple
nohup perl auto_eggNOG.pl -workdir ./ -eggNOG -GO -KEGG -COG

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/../../lib";
use GACP qw(parse_config);
use Cwd qw(getcwd abs_path);

my ($workdir,$eggNOG,$GO,$KEGG,$COG,$help);
GetOptions(
	"workdir:s"=>\$workdir,
	"eggNOG"=>\$eggNOG,
    "GO"=>\$GO,
    "KEGG"=>\$KEGG,
    "COG"=>\$COG,
	"Help"=>\$help,
);
die `pod2text $0` if(@ARGV==0 || $help);

my $fasta = shift;
$workdir ||= "./";
$fasta = abs_path($fasta);
$workdir = abs_path($workdir);
my $filename = basename($fasta);

my $config_file = "$Bin/../../config.txt";
my $func_bin = parse_config($config_file,"func_bin");
my $go_bin = "$Bin/go";
my $kegg_bin = "$Bin/kegg";
my $eggNOG_path = parse_config($config_file,"eggNOG");
my $GO_path  = parse_config($config_file,"GO");
my $KEGG_path  = parse_config($config_file,"KEGG");
my $COG_path = parse_config($config_file,"COG");
my $activate = parse_config($config_file,"activate");
my $convert = parse_config($config_file,"convert");


if (defined $eggNOG){
	&mkdir_chdir("$workdir/eggNOG");
	print "Start eggNOG and COG function annotation \nMake directory  $workdir/eggNOG && cd  $workdir/eggNOG ..\n";

	open SH_eggNOG, ">$workdir/eggNOG/$filename.eggNOG.sh" || die $!;
	print SH_eggNOG "#!/bin/bash\n";
	print SH_eggNOG "## eggNOG analysis ###\n";
	print SH_eggNOG "source $activate denovo\n";
	print SH_eggNOG "emapper.py -m diamond -i $fasta --output $workdir/eggNOG/out --override --cpu 10 --data_dir $eggNOG_path \n";
	print SH_eggNOG "source $activate base_env \n";
	print SH_eggNOG "python $Bin/EggnogParser_change.py -f $fasta -n $workdir/eggNOG/out.emapper.annotations -o $workdir/eggNOG \n";
	print SH_eggNOG "Rscript $Bin/draw_COG.R $workdir/COG/all.COG.class.txt $workdir/COG \n";
	close SH_eggNOG;
	system("bash $workdir/eggNOG/$filename.eggNOG.sh 1>$workdir/eggNOG/$filename.eggNOG.sh.e 2>$workdir/eggNOG/$filename.eggNOG.sh.o") == 0 || die $!;
}

if (defined $GO){
	&mkdir_chdir("$workdir/GO_eggNOG");
	print "Start GO function annotation \nMake directory  $workdir/GO_eggNOG && cd  $workdir/GO_eggNOG ..\n";

	open SH_GO, ">$workdir/GO_eggNOG/$filename.GO.sh" || die $!;
	print SH_GO "#!/bin/bash\n";
	print SH_GO "### GO analysis ###\n";
	print SH_GO "cut -f 1,10 $workdir/eggNOG/out.emapper.annotations | perl $go_bin/get_go_anno.pl - > $workdir/GO_eggNOG/eggNOG2GO.anno && \\\n";
	print SH_GO "mkdir -p $workdir/GO_eggNOG/data && \\\n";
	print SH_GO "perl $go_bin/annot2goa.pl $GO_path/go-basic.obo $workdir/GO_eggNOG/eggNOG2GO.anno $workdir/GO_eggNOG/data/species && \\\n";
	print SH_GO "cat $workdir/GO_eggNOG/data/species.[CFP] | cut -f 2 | sort -u > $workdir/GO_eggNOG/gene.list && \\\n";
	print SH_GO "perl $go_bin/drawGO.pl -list $workdir/GO_eggNOG/gene.list -goprefix $workdir/GO_eggNOG/data/species -goclass $GO_path/go.class -outprefix $filename \n";
	close SH_GO;
	system("bash $workdir/GO_eggNOG/$filename.GO.sh 1>$workdir/GO_eggNOG/$filename.GO.sh.e 2>$workdir/GO_eggNOG/$filename.GO.sh.o") == 0 || die $!;
}

if (defined $KEGG){
	&mkdir_chdir("$workdir/KEGG");
	print "Start KEGG function annotation \nMake directory  $workdir/KEGG && cd  $workdir/KEGG ..\n";

	open SH_KEGG, ">$workdir/KEGG/$filename.KEGG.sh" || die $!;
	print SH_KEGG "#!/bin/bash\n";
	print SH_KEGG "### KEGG analysis ###\n";
	print SH_KEGG "cut -f 1,12 $workdir/eggNOG/out.emapper.annotations | perl $kegg_bin/get_kegg_desc.pl $KEGG_path/ko.txt.desc - > $workdir/KEGG/$filename.kegg.anno \n";
	print SH_KEGG "awk '\$2!=\"-\"' $workdir/KEGG/$filename.kegg.anno > $workdir/KEGG/$filename.kegg.txt\n";
	print SH_KEGG "cat  $workdir/KEGG/$filename.kegg.txt |cut -f 2 |sed 's/,/\\n/g' |grep -v '^\$' |grep -v '-' | sed 's/ko://g' | perl $Bin/fishInWinter.pl - $KEGG_path/ko.txt.level |cut -f 2,3 > $workdir/KEGG/R.input \n";
	print SH_KEGG "if [ -s $workdir/KEGG/R.input ];then \n";
    print SH_KEGG "\t Rscript $Bin/draw_kegg.R -i $workdir/KEGG/R.input -o $workdir/KEGG/$filename.KEGG.pdf \n";
    print SH_KEGG "\t $convert -density 300 $workdir/KEGG/$filename.KEGG.pdf $workdir/KEGG/$filename.KEGG.png \n";
    print SH_KEGG "else\n\techo \"No KEGG anno\"\n";
    print SH_KEGG "fi\n";
	
	close SH_KEGG;
	system("bash $workdir/KEGG/$filename.KEGG.sh 1>$workdir/KEGG/$filename.KEGG.sh.e 2>$workdir/KEGG/$filename.KEGG.sh.o ") == 0 || die $!;
}

# if (defined $COG){
# 	&mkdir_chdir("$workdir/COG");
# 	print "Start COG function annotation \nMake directory  $workdir/COG && cd  $workdir/COG ..\n";

# 	open SH_COG, ">$workdir/COG/$filename.COG.sh" || die $!;
# 	print SH_COG "#!/bin/bash\n";
# 	print SH_COG "## COG analysis ###\n";
# 	print SH_COG "cut -f 1,7 $workdir/eggNOG/out.emapper.annotations | perl $func_bin/get_cog_desc.pl $COG_path/cog.db.class.desc - $workdir/COG/$filename.cog.xls $workdir/COG/cog.num.stat && \\\n";
# 	print SH_COG "perl $func_bin/cog_R.pl -catalog $workdir/COG/cog.num.stat -outdir ./ -sample $filename -Rscript Rscript -convert convert\n";
# 	close SH_COG;
# 	system("bash $workdir/COG/$filename.COG.sh 1>$workdir/COG/$filename.COG.sh.e 2>$workdir/COG/$filename.COG.sh.o") == 0 || die $!;
# }

#####################################sub#################################
sub mkdir_chdir{
	my $dir_name = shift;
	system("rm -rf $dir_name") == 0 || die $! if (-d $dir_name);
	system("mkdir $dir_name") == 0 || die $!;
	chdir $dir_name or die ("ERROR: Cannot change directory to $dir_name: $!\n");
}
