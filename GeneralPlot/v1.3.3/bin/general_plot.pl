#!/Bio/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
# Program:       general_plot.pl
# Author:        Liu yubin
# Date:          Thu 09 Apr 2020 04:04:55 PM CST

use lib "$Bin/../lib";
use GeneralPlot;

general_plot(@ARGV);

