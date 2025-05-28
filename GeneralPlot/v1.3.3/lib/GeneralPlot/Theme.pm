package GeneralPlot::Theme;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.13 17:53:51          |
#--------------------------------------------------#
use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Config::General;
use parent 'Exporter';
our @EXPORT = (
    'sbv_venn_theme_default',      'ggplot2_line_theme_default', 'ggplot2_pie_theme_default', 
    'ggplot2_bar_theme_default',   'ggplot2_hist_theme_default', 'ggplot2_heatmap_theme_default',
    'ggplot2_heatmap_theme_ggcor', 'ggplot2_dot_theme_default',  'ggplot2_dot_theme_manhattan',
    'ggplot2_box_theme_default',
);
our @EXPORT_OK = qw();

use lib "$Bin";
use lib "$Bin/..";
use GeneralPlot::Debug;
use GeneralPlot::Envs;

sub sbv_venn_theme_default {

    my %confs = @_;

    $confs{'-circle_stroke'} ||= '\#000000';
    $confs{'-circle_stroke-width'} ||= 0;
    $confs{'-circle_fill-opacity'} ||= 0.8;

    $confs{'-ellipse_stroke'} ||= '\#000000';
    $confs{'-ellipse_stroke-width'} ||= 0;
    $confs{'-ellipse_fill-opacity'} ||= 0.8;

    $confs{'-path_stroke'} ||= '\#000000';
    $confs{'-path_stroke-width'} ||= 0;
    $confs{'-path_fill-opacity'} ||= 0.8;

    $confs{'-text_fill'} ||= '\#000000';
    $confs{'-text_font-family'} ||= 'arial';
    $confs{'-text_font-style'} ||= 'normal';
    $confs{'-text_font-weight'} ||= 'bold';
    $confs{'-text_font-size'} ||= 22;
    $confs{'-text_stroke-width'} ||= 0;
    
    $confs{'-text.title_fill'} ||= '\#000000';
    $confs{'-text.title_font-family'} ||= 'arial';
    $confs{'-text.title_font-style'} ||= 'normal';
    $confs{'-text.title_font-weight'} ||= 'bold';
    $confs{'-text.title_font-size'} ||= 24;
    $confs{'-text.title_stroke-width'} ||= 0;

    $confs{'-text.label_fill'} ||= '\#000000';
    $confs{'-text.label_font-family'} ||= 'arial';
    $confs{'-text.label_font-style'} ||= 'normal';
    $confs{'-text.label_font-weight'} ||= 'bold';
    $confs{'-text.label_font-size'} ||= 24;
    $confs{'-text.label_stroke-width'} ||= 0;

    my $conf = <<CONF;
<styles>
    <circle>
        stroke = $confs{'-circle_stroke'}
        stroke-width = $confs{'-circle_stroke-width'}
        fill-opacity = $confs{'-circle_fill-opacity'}
    </circle>

    <ellipse>
        stroke = $confs{'-ellipse_stroke'}
        stroke-width = $confs{'-ellipse_stroke-width'}
        fill-opacity = $confs{'-ellipse_fill-opacity'}
    </ellipse>

    <path>
        stroke = $confs{'-path_stroke'}
        stroke-width = $confs{'-path_stroke-width'}
        fill-opacity = $confs{'-path_fill-opacity'}
    </path>

    <text>
        fill = $confs{'-text_fill'}
        font-family = $confs{'-text_font-family'}
        font-style = $confs{'-text_font-style'}
        font-weight = $confs{'-text_font-weight'}
        font-size = $confs{'-text_font-size'}
        stroke-width = $confs{'-text_stroke-width'}
    </text>

    <text.title>
        fill = $confs{'-text.title_fill'}
        font-family = $confs{'-text.title_font-family'}
        font-style = $confs{'-text.title_font-style'}
        font-weight = $confs{'-text.title_font-weight'}
        font-size = $confs{'-text.title_font-size'}
        stroke-width = $confs{'-text.title_stroke-width'}
    </text.title>

    <text.label>
        fill = $confs{'-text.label_fill'}
        font-family = $confs{'-text.label_font-family'}
        font-style = $confs{'-text.label_font-style'}
        font-weight = $confs{'-text.label_font-weight'}
        font-size = $confs{'-text.label_font-size'}
        stroke-width = $confs{'-text.label_stroke-width'}
    </text.label>
</styles>
CONF

    return $conf;

}

sub ggplot2_line_theme_default {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- line_theme_default()
RSC

    return $rscript;

}

sub ggplot2_pie_theme_default {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- pie_theme_default()
RSC

    return $rscript;

}

sub ggplot2_bar_theme_default {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- bar_theme_default()
RSC

    return $rscript;

}

sub ggplot2_hist_theme_default {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- hist_theme_default()
RSC

    return $rscript;

}

sub ggplot2_heatmap_theme_default {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- heatmap_theme_default()
RSC

    return $rscript;

}

sub ggplot2_heatmap_theme_ggcor {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- heatmap_theme_ggcor()
RSC

    return $rscript;

}

sub ggplot2_dot_theme_default {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- dot_theme_default()
RSC

    return $rscript;

}


sub ggplot2_dot_theme_manhattan {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- dot_theme_manhattan()
RSC

    return $rscript;

}

sub ggplot2_box_theme_default {

    my $rscript = <<RSC;
source("$SOFT{all_theme}", chdir = T)
mytheme <- box_theme_default()
RSC

    return $rscript;

}

return 1;

