package GeneralPlot::Image::Line;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.10 16:10:47          |
#--------------------------------------------------#
use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Config::General;
use parent 'Exporter';
our @EXPORT = qw();
our @EXPORT_OK = qw();

use lib "$RealBin";
use lib "$RealBin/..";
use lib "$RealBin/../..";

use GeneralPlot::Envs;
use GeneralPlot::Debug;
use GeneralPlot::Data;
use GeneralPlot::Theme;
use GeneralPlot::Colors;

#-------------------------------------------------------------------------------
#  new a line object 
#-------------------------------------------------------------------------------
sub new {

    my ($class, %opts) = @_;

    my $self = {
        'type'   => 'line',
        'method' => 'none',
        'obj'    => 0,
        'save'   => 0,
        'run'    => 0,
    };

    bless $self, $class;

    return $self;

}

#-------------------------------------------------------------------------------
#  general a line chart
#-------------------------------------------------------------------------------
sub plot {

    my ($self, %opts) = @_;

    $opts{'-method'} ||= "line";

    if(defined $opts{'-method'} and $opts{'-method'} eq 'line') {
        $self->{method} = $opts{'-method'};
        plot_line($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'cumulative_line') {
        $self->{method} = $opts{'-method'};
        plot_cumulative($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'stauration') {
        $self->{method} = $opts{'-method'};
        plot_stauration($self, %opts);
    }
    else {
        ERROR("$self->{type}", "check the -method [$opts{'-method'}]. ");
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general line chart from a format table directly
#-------------------------------------------------------------------------------
sub plot_line {

    my ($self, %opts) = @_;

    if(defined $opts{'-outprefix'}) {
        $opts{'-outdir'}  ||= dirname($opts{'-outprefix'});
        $opts{'-outname'} ||= basename($opts{'-outprefix'});
    }
    !defined $self->{pic}->{file} and !defined $opts{'-file'} and 
    ERROR("$self->{type}", "<-file> is not defined. ");
    !defined $self->{pic}->{outdir} and !defined $opts{'-outdir'} and 
    ERROR("$self->{type}", "<-outdir> is not defined. ");
    !defined $self->{pic}->{outname} and !defined $opts{'-outname'} and
    ERROR("$self->{type}", "<-outname> is not defined. ");

    defined $opts{'-file'}    and $self->{pic}->{file}    = rel2abs($opts{'-file'});
    defined $opts{'-outdir'}  and $self->{pic}->{outdir}  = rel2abs($opts{'-outdir'});
    defined $opts{'-outname'} and $self->{pic}->{outname} = $opts{'-outname'};

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";

    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{title}   = $opts{'-title'}   ? $opts{'-title'}   : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}    ? $opts{'-xlab'}    : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}    ? $opts{'-ylab'}    : "";
    $self->{pic}->{legend}  = $opts{'-legend'}  ? $opts{'-legend'}  : "T";
    $self->{pic}->{width}   = $opts{'-width'}   ? $opts{'-width'}   : 7;
    $self->{pic}->{height}  = $opts{'-height'}  ? $opts{'-height'}  : 6;
    $self->{pic}->{omitNA}  = $opts{'-omitNA'}  ? $opts{'-omitNA'}  : "T";

    # get colors
    my $color = "";
    my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix");
    my $sample_cnt = $ncol - 1;
    if($sample_cnt < 1) {
        ERROR("line", "no sample found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("line", "default", $sample_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $sample_cnt) {
            $color = join(",", map {"'$_'"} @$colors[0 .. ($sample_cnt - 1)]);
        }
        else {
            $color = join(",", map { "'$_'" } @$colors);
            $color = "colorRampPalette(c($color))($sample_cnt)";
        }
    }

    if($self->{pic}->{colnpg} =~ /^T$|^True$/i) {
        my $colors = fetch_color(-type => "ggsci", -tag => 'npg', -num => 0);
        $color = join(",", map {"'$_'"} @$colors);
    }
    elsif($self->{pic}->{collanc} =~ /^T$|^True$/i) {
        my $colors = fetch_color(-type => "ggsci", -tag => 'lancet', -num => 0);
        $color = join(",", map {"'$_'"} @$colors);
    }
    elsif($self->{pic}->{collist} ne 'none') {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            if(scalar(@colors) != $sample_cnt) {
                my $col_cnt = scalar(@colors);
                warn("Insufficient colors. $sample_cnt needed but only $col_cnt provided. \n");
            }
            $color = join(",", map {"'$_'"} @colors);
        }
        elsif($sample_cnt == 1) {
            $color = "'$self->{pic}->{collist}'";
        }
    }

    # width adjust 
    if($sample_cnt <= 20) {
        $self->{pic}->{width} = $opts{'-width'} ? $opts{'-width'} : 7;
    }
    elsif($sample_cnt <= 40) {
        $self->{pic}->{width} = $opts{'-width'} ? $opts{'-width'} : 8;
    }
    elsif($sample_cnt <= 60) {
        $self->{pic}->{width} = $opts{'-width'} ? $opts{'-width'} : 9;
    }
    else {
        $self->{pic}->{width} = $opts{'-width'} ? $opts{'-width'} : 7 + $sample_cnt / 20;
    }
    if($self->{pic}->{legend} !~ /^T$|^True$/i) {
        $self->{pic}->{width} = $opts{'-width'} ? $opts{'-width'} : 7;
    }

    # get theme
    my $header = ggplot2_head();
    my $theme  = ggplot2_line_theme_default();

    my $rscript = <<RSC;
$header
$theme
source("$SOFT{ggplot2}")

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

title  <- "$self->{pic}->{title}"
xlab   <- "$self->{pic}->{xlab}"
ylab   <- "$self->{pic}->{ylab}"
islegend <- as.logical("$self->{pic}->{legend}")
isomitNA <- as.logical("$self->{pic}->{omitNA}")

width  <- $self->{pic}->{width}
height <- $self->{pic}->{height}
colors <- c($color)

data <- read.table(file, header = $self->{pic}->{header}, sep = "\\t", quote = "", check.names = F)
colnames(data)[1] = "x__index"
data.melt <- melt(data, id.vars = c("x__index"))
sample.num <- length(unique(data.melt\$variable))
data.melt\$variable <- factor(
    data.melt\$variable, levels = unique(data.melt\$variable), ordered = T)
if(is.numeric(data\$x__index) == FALSE) {
    data.melt\$x__index <- factor(
        data.melt\$x__index, levels = unique(data.melt\$x__index), ordered = T)
}

if(isomitNA) {
    data.melt <- na.omit(data.melt)
}

p <- ggplot(data.melt, aes(x = x__index, y = value, color = variable, group = variable)) +
geom_line(size = 1, alpha = 0.9)

# we assume that the data of y-axis is continuous data
if(is.numeric(data\$x__index) == TRUE) {
    x_min <- min(data\$x__index)
    x_max <- max(data\$x__index)
    y_min <- min(data.melt\$value)
    y_max <- max(data.melt\$value) * 1.12
    if(x_min >= 0) {
        x_min = 0
    }
    if(y_min >= 0) {
        y_min = 0
    }

    p <- p + 
    scale_x_continuous(expand = c(0.02, 0.02), breaks = pretty_breaks(5)) +
    scale_y_continuous(expand = c(0.01, 0.01), breaks = pretty_breaks(5), limits = c(y_min, y_max))
} else {
    y_min <- min(data.melt\$value)
    y_max <- max(data.melt\$value) * 1.12
    if(y_min >= 0) {
        y_min = 0
    }

    width = 8
    x_idx.num <- length(unique(data.melt\$x__index))
    x_expand = c(0.03, 0.03)
    if(x_idx.num <= 10) {
        width = width
    } else if(x_idx.num <= 30) {
        width = width + 2
        x_expand = c(0.02, 0.02)
    } else if(x_idx.num <= 100) {
        width = width + 2 + 8 * (x_idx.num - 30) / 70
        x_expand = c(0.015, 0.015)
    } else {
        width = 20
        x_expand = c(0.012, 0.012)
    }
    axis.size = AxisTextSize(data.melt\$x__index, width)

    p <- p + scale_x_discrete(expand = x_expand) +
    scale_y_continuous(expand = c(0.01, 0.01), breaks = pretty_breaks(5), limits = c(y_min, y_max))
}

if(length(colors) >= sample.num) {
    colors <- colors[1:sample.num]
} else {
    colors <- colorRampPalette(colors)(sample.num)
}
p <- p + scale_color_manual(values = colors) +
labs(x = xlab, y = ylab, title = title) + 
mytheme

if(is.numeric(data\$x__index) == FALSE) {
    p <- p + theme(axis.text.x = element_text(size = axis.size)) + 
    xlab_angle(unique(data\$x), width = width)
}
if(sample.num == 1) {
    width = 6
    p <- p + theme(legend.position = "none")
}
if(islegend != T) {
    width = 6
    p <- p + theme(legend.position = "none")
}

ggsave(filename = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "line";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
# plot the cumulative lines with a matrix
#-------------------------------------------------------------------------------
sub plot_cumulative {

    my ($self, %opts) = @_;
    
    if(defined $opts{'-outprefix'}) {
        $opts{'-outdir'}  ||= dirname($opts{'-outprefix'});
        $opts{'-outname'} ||= basename($opts{'-outprefix'});
    }

    !defined $self->{pic}->{file} and !defined $opts{'-file'} and
    ERROR("$self->{type}", "<-file> is not defined. ");
    !defined $self->{pic}->{outdir} and !defined $opts{'-outdir'} and
    ERROR("$self->{type}", "<-outdir> is not defined. ");
    !defined $self->{pic}->{outname} and !defined $opts{'-outname'} and
    ERROR("$self->{type}", "<-outname> is not defined. ");

    defined $opts{'-file'}    and $self->{pic}->{file}    = rel2abs($opts{'-file'});
    defined $opts{'-outdir'}  and $self->{pic}->{outdir}  = rel2abs($opts{'-outdir'});
    defined $opts{'-outname'} and $self->{pic}->{outname} = $opts{'-outname'};

    $self->{pic}->{group}   = $opts{'-group'}   ? rel2abs($opts{'-group'}) : "none";
    $self->{pic}->{scale}   = $opts{'-scale'}   ? $opts{'-scale'} : "none";
    $self->{pic}->{filter0} = $opts{'-filter0'} ? $opts{'-filter0'} : "FALSE";
    $self->{pic}->{namefix} = $opts{'-namefix'} ? $opts{'-namefix'} : "none";
    $self->{pic}->{color}   = $opts{'-color'}   ? $opts{'-color'} : "none";
    $self->{pic}->{title}   = $opts{'-title'}   ? $opts{'-title'} : "none";
    $self->{pic}->{x_lab}   = $opts{'-x_lab'}   ? $opts{'-x_lab'} : "none";
    $self->{pic}->{y_lab}   = $opts{'-y_lab'}   ? $opts{'-y_lab'} : "none";
    $self->{pic}->{size}    = $opts{'-size'}    ? $opts{'-size'} : 1;
    $self->{pic}->{alpha}   = $opts{'-alpha'}   ? $opts{'-alpha'} : 0.9;
    $self->{pic}->{width}   = $opts{'-width'}   ? $opts{'-width'} : 8;
    $self->{pic}->{height}  = $opts{'-height'}  ? $opts{'-height'} : 7;
    $self->{pic}->{out_format} = $opts{'-out_format'} ? $opts{'-out_format'} : "pic_data";

    $self->{pic}->{scale}   = "'$self->{pic}->{scale}'";
    $self->{pic}->{namefix} = "'$self->{pic}->{namefix}'";
    $self->{pic}->{color}   = "'$self->{pic}->{color}'";
    $self->{pic}->{title}   = "'$self->{pic}->{title}'";
    $self->{pic}->{x_lab}   = "'$self->{pic}->{x_lab}'";
    $self->{pic}->{y_lab}   = "'$self->{pic}->{y_lab}'";

    my $shell = "export LANG='zh_CN.UTF-8' ; mkdir -p $self->{pic}->{outdir} ; ".
    "$SOFT{Rscript} $SOFT{cumulative_line} $self->{pic}->{file} ".
    "$self->{pic}->{outdir}/$self->{pic}->{outname} $self->{pic}->{group} ".
    "$self->{pic}->{filter0} $self->{pic}->{scale} $self->{pic}->{namefix} ".
    "$self->{pic}->{color} $self->{pic}->{title} $self->{pic}->{x_lab} ".
    "$self->{pic}->{y_lab} $self->{pic}->{size} $self->{pic}->{alpha} ".
    "$self->{pic}->{width} $self->{pic}->{height} $self->{pic}->{out_format}";
    
    $self->{rscript} = "";
    $self->{command} = $shell;
    $self->{method} = "cumulative_line";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  specific chart process for gene stauration
#-------------------------------------------------------------------------------
#/Bio/Bin/pipeline//System_Programs/Small_pipe/v1.3/R_plot/src/stauration_analysis.R
sub plot_stauration {

    my ($self, %opts) = @_;


    return $self;

}

return 1;

