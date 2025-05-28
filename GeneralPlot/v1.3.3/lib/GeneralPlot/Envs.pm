package GeneralPlot::Envs;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.05.06 18:03:47          |
#--------------------------------------------------#
use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Config::General;
use base 'Exporter';
our @EXPORT = qw($APP_NAME $VERSION %SOFT ggplot2_head load_opts save run);
our @EXPORT_OK = qw();

use lib "$RealBin";
use lib "$RealBin/..";

use GeneralPlot::Debug;

my $bin_path = dirname(rel2abs(__FILE__));
my $src_path = $bin_path . "/../../src";
$src_path = abs_path($src_path);

#-------------------------------------------------------------------------------
# Global Variable
#-------------------------------------------------------------------------------
our $APP_NAME = "GeneralPlot";
our $AUTHOR   = "genedenovo";
our $VERSION  = "v1.3.0";

our %SOFT = (

    'sh'              => 'sh',
    'bash'            => 'bash',
    'qsub'            => 'qsub',
    
    'perl'            => '/usr/bin/perl',
    'Rscript'         => '/data1/NFS/home/rcb/software/Miniconda/miniconda3/envs/R-4.1.1/bin/Rscript',
    'convert'         => 'convert',
    
    'sbv'             => "$src_path/SBV/bin/sbv.pl",
    'colors_run'      => "$src_path/Rlib/Colors/Colors_run.r",
    'get_color_pal'   => "$src_path/Rlib/get_color_Palette.r",
    'all_theme'       => "$src_path/Rlib/all.theme.r",
    'ggplot2'         => "$src_path/Rlib/ggplot2.r",
    'split_violin'    => "$src_path/Rlib/split_violin/split_violin_ggplot.r",
    'scatter_3d'      => "$src_path/Rlib/scatter_3d/scatter3D_plot.r",
    'polar_bar'       => "$src_path/Rlib/polar_bar/polar_bar.r", 
    'sankey_plot'     => "$src_path/Rlib/sankey_plot/sankey_plot.R",
    'scatter_bubble'  => "$src_path/Rlib/scatter_bubble/scatter_bubble.R",
    'lollipop'        => "$src_path/Rlib/lollipop/lollipop.R",
    'cumulative_line' => "$src_path/Rlib/cumulative_line/cumulative_line.R",

);

# %SOFT = change_env_soft(%SOFT);

#-------------------------------------------------------------------------------
# sub functions
#-------------------------------------------------------------------------------
sub change_env_soft {

    my %SOFT = @_;

    if($ENV{HOSTNAME} ne "genedenovo.org") {
        if($ENV{HOSTNAME} !~ /compute-/) {
            $SOFT{'Rscript'} = "/public/home/pangrui/RCB/softwares/miniconda3/miniconda3-4.7.12/envs/R-4.1.1/bin/Rscript";
        }
    }

    return %SOFT;
}

sub ggplot2_head {

    my $rscript = <<RSC;
#!$SOFT{'Rscript'}
library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
library(RColorBrewer)
#library(showtext)
RSC

    return $rscript;

}

#-------------------------------------------------------------------------------
#  load the opts
#-------------------------------------------------------------------------------
sub load_opts {

    my ($self, %opts) = @_;

    for my $opt (keys %opts) {
        $self->{$opt} = $opts{$opt};
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general scripts
#-------------------------------------------------------------------------------
sub save {

    my ($self, %opts) = @_;

    $opts{'-clean'} ||= "F";    # T | F
    $self->{convert_density} = defined $self->{convert_density} ? $self->{convert_density} : 
                               defined $opts{'-convert_density'} ? $opts{'-convert_density'} : 300;

    if($self->{obj} == 0) {
        ERROR("save", "you should run $self->{type}() first. ");
    }
    `mkdir -p $self->{pic}->{'outdir'}` unless(-d $self->{pic}->{'outdir'});

    if(exists $self->{rscript} and $self->{rscript} ne "") {
        open RSC, "> $self->{pic}->{outdir}/$self->{pic}->{outname}.r" or ERROR("save", $!);
        print RSC $self->{rscript};
        close RSC;
        $self->{tmp_file}->{rscript} = "$self->{pic}->{outdir}/$self->{pic}->{outname}.r";
    }
    elsif(exists $self->{conf}) {
        open CFG, "> $self->{pic}->{outdir}/$self->{pic}->{outname}.conf" or ERROR("save", "$!");
        print CFG "$self->{conf}";
        close CFG;
        $self->{tmp_file}->{conf} = "$self->{pic}->{outdir}/$self->{pic}->{outname}.conf";
    }

    open SH, "> $self->{pic}->{outdir}/$self->{pic}->{outname}.sh" or ERROR("save", "$!");
    print SH "#!$SOFT{bash} \n";
    print SH "set -e \n";
    if(exists $self->{rscript} and $self->{rscript} ne "") {
        if($self->{type} eq 'heatmap' and $self->{method} eq 'cluster') {
            print SH "$SOFT{'Rscript'} $self->{tmp_file}->{rscript} \n";
        }
        elsif($self->{type} eq 'dot' and $self->{method} eq 'manhattan') {
            print SH "$SOFT{'Rscript'} $self->{tmp_file}->{rscript} \n";
        }
        else {
            print SH "$SOFT{Rscript} $self->{tmp_file}->{rscript} \n";
            print SH "$SOFT{convert} -density $self->{convert_density} $self->{pic}->{outdir}/".
            "$self->{pic}->{outname}.pdf $self->{pic}->{outdir}/$self->{pic}->{outname}.png \n";
        }
    }
    elsif(exists $self->{conf}) {
        print SH "$SOFT{perl} $SOFT{sbv} venn -conf $self->{tmp_file}->{conf} \n";
        print SH "$SOFT{convert} -density $self->{convert_density} $self->{pic}->{outdir}/".
        "$self->{pic}->{outname}.svg $self->{pic}->{outdir}/$self->{pic}->{outname}.png \n";
    }
    elsif(exists $self->{command}) {
        print SH "$self->{command}\n";
        #print SH "$SOFT{convert} -density $self->{convert_density} $self->{pic}->{outdir}/".
        #"$self->{pic}->{outname}.pdf $self->{pic}->{outdir}/$self->{pic}->{outname}.png \n";
        print SH "if [ -s $self->{pic}->{outdir}/$self->{pic}->{outname}.pdf ] ; then ".
        "$SOFT{convert} -density $self->{convert_density} $self->{pic}->{outdir}/".
        "$self->{pic}->{outname}.pdf $self->{pic}->{outdir}/$self->{pic}->{outname}.png ; fi \n";
    }
    $self->{sh_file} = "$self->{pic}->{outdir}/$self->{pic}->{outname}.sh";
    $self->{tmp_file}->{shell} = "$self->{pic}->{outdir}/$self->{pic}->{outname}.sh";

    if($opts{'-clean'} =~ /^T$|^True$/i) {
        if(defined $self->{tmp_file}) {
            for my $key (sort keys %{$self->{tmp_file}}) {
                print SH "rm $self->{tmp_file}->{$key} \n";
            }
        }
        if(defined $self->{rm_file}) {
            for my $key (sort keys %{$self->{rm_file}}) {
                print SH "rm $self->{rm_file}->{$key} \n";
            }
        }
    }
    close SH;

    $self->{save} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  run scripts
#-------------------------------------------------------------------------------
sub run {

    my ($self, %opts) = @_;

    timeLOG("GeneralPlot [$self->{type}::$self->{method}] '$self->{pic}->{outname}' start ... ");

    if($self->{save} == 0) {
        $self->save(%opts);
    }

    system("$SOFT{sh} $self->{sh_file}");
    $self->{run} = 1;

    timeLOG("GeneralPlot [$self->{type}::$self->{method}] '$self->{pic}->{outname}' finished. ");

    return 0;

}

return 1;

