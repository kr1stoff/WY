package GeneralPlot;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.09 16:06:12          |
#--------------------------------------------------#
use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Config::General;
use parent 'Exporter';
our @EXPORT = qw(general_plot);
our @EXPORT_OK = qw();

use lib "$Bin";

use GeneralPlot::Envs;
use GeneralPlot::Debug;
use GeneralPlot::Data;
use GeneralPlot::Stat;
use GeneralPlot::Image;
use GeneralPlot::Theme;
use GeneralPlot::Colors;
use GeneralPlot::Help;

#-------------------------------------------------------------------------------
#  main function
#-------------------------------------------------------------------------------

sub general_plot {

    my $type = shift;

    ## a rough argv check without using GetOptions
    if(@_ == 0) {
        usage();
        exit 1;
    }
    else {
        for my $i (@_) {
            if($i =~ /^-help$|^-h$/) {
                my $flag = type_check($type);
                if($flag == 0) {
                    sub_usage($type);
                }
                elsif($flag == 1) {
                    usage();
                }
                exit 0;
            }
            elsif($i =~ /^-version$|^-v$/) {
                print "** Version: $APP_NAME: $VERSION ** \n";
                exit 0;
            }
        }
    }

    ## recheck type
    my $flag = type_check($type);
    if($flag != 0) {
        usage();
        ERROR('type_check', "[$type] is not vaild type. ");
    }

    ## plot 
    my %opts = @_;
    my $p = GeneralPlot::Image->new();
    $p->plot(-graph => $type, %opts);

}

return 1;

