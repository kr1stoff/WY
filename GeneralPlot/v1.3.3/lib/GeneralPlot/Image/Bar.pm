package GeneralPlot::Image::Bar;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.10 16:14:24          |
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
use GeneralPlot::Colors;
use GeneralPlot::Theme;
use GeneralPlot::Data;

#-------------------------------------------------------------------------------
#  new a  object 
#-------------------------------------------------------------------------------
sub new {

    my ($class, %opts) = @_;

    my $self = {
        'type'   => 'bar',
        'method' => 'none',
        'obj'    => 0,
        'save'   => 0,
        'run'    => 0,
    };

    bless $self, $class;

    return $self;

}

sub plot {

    my ($self, %opts) = @_;

    $opts{'-method'} ||= "bar";

    if(defined $opts{'-method'} and $opts{'-method'} eq 'bar') {
        $self->{method} = $opts{'-method'};
        bar($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'errorbar') {
        $self->{method} = $opts{'-method'};
        error_bar($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'polarbar') {
        $self->{method} = $opts{'-method'};
        polar_bar($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'sankey') {
        $self->{method} = $opts{'-method'};
        sankey($self, %opts);
    }
    else {
        ERROR("$self->{type}", "check the -method [$opts{'-method'}]. ");
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general a single bar chart from a format table directly
#-------------------------------------------------------------------------------
sub bar {

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

    $self->{pic}->{header}   = $opts{'-header'}   ? $opts{'-header'}   : "T";
    $self->{pic}->{bracket}  = $opts{'-bracket'}  ? $opts{'-bracket'}  : "T";

    $self->{pic}->{barfill}   = $opts{'-barfill'}   ? $opts{'-barfill'}   : "F";
    $self->{pic}->{bardodge}  = $opts{'-bardodge'}  ? $opts{'-bardodge'}  : "F";
    $self->{pic}->{mulcolor}  = $opts{'-mulcolor'}  ? $opts{'-mulcolor'}  : "F";
    $self->{pic}->{coordflip} = $opts{'-coordflip'} ? $opts{'-coordflip'} : "F";

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";    # 'red,blue'

    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{label}   = $opts{'-label'}  ? $opts{'-label'}  : "F";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 6;

    # color
    my $color = "";
    my $label_cnt = 0;
    my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix");
    my ($sample_cnt, $attr_cnt) = ($nrow, $ncol - 1);
    $sample_cnt-- if($self->{pic}->{header} =~ /^T$|^True$/i);
    
    if($self->{pic}->{mulcolor} =~ /^T$|^True$/i and $attr_cnt == 1) {
        $label_cnt = $sample_cnt;
    }
    else {
        $label_cnt = $attr_cnt;
    }

    if($label_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("bar", "default", $label_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if($self->{pic}->{bardodge} !~ /^T$|^True$/i and $self->{pic}->{colset} eq 'none') {
            #if($label_cnt > 10) {
            #    $colors = fetch_color(-type => "bar", -tag => 'stack', -num => 0);
            #}
            #elsif($label_cnt > 1) {
                $colors = fetch_color(-type => "bar", -tag => 'stack', -num => $num);
            #}
        }
        if(scalar(@$colors) >= $label_cnt) {
            $color = join(",", map { "'$_'" } @$colors[0 .. ($label_cnt - 1)]);
        }
        else {
            $color = join(",", map { "'$_'" } @$colors);
            $color = "colorRampPalette(c($color))($label_cnt)";
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
    elsif($self->{pic}->{collist} ne "none") {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            $color = join(",", map {"'$_'"} @colors);
        }
        else {
            $color = "'$self->{pic}->{collist}'";
        }
    }
    
    # default label settings 
    #if($sample_cnt <= 10 and $attr_cnt == 1) {
    #    $self->{pic}->{label} = $opts{'-label'} ? $opts{'-label'} : "T";
    #}
    #elsif($self->{pic}->{bardodge} =~ /^T$|^True$/i and $attr_cnt <= 2 and $sample_cnt <= 10) {
    #    $self->{pic}->{label} = $opts{'-label'} ? $opts{'-label'} : "T";
    #}

    # theme
    my $header = ggplot2_head();
    my $theme  = ggplot2_bar_theme_default();

    my $rscript = <<RSC;
$header
$theme
source("$SOFT{'ggplot2'}")

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

title   <- "$self->{pic}->{title}"
xlab    <- "$self->{pic}->{xlab}"
ylab    <- "$self->{pic}->{ylab}"
islabel   <- as.logical("$self->{pic}->{label}")
islegend  <- as.logical("$self->{pic}->{legend}")
isbracket <- as.logical("$self->{pic}->{bracket}")
iscoordflip <- as.logical("$self->{pic}->{coordflip}")

width   <- $self->{pic}->{width}
height  <- $self->{pic}->{height}
colors <- c($color)

data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
colnames(data)[1] = "x_index"
#is.x.num = T
#if(is.numeric(data\$x_index) == F) {
#    data\$x_index <- factor(data\$x_index, levels = unique(data\$x_index), ordered = T)
#    is.x.num = F
#}
data\$x_index <- factor(data\$x_index, levels = unique(data\$x_index), ordered = T)
data.melt <- melt(data, id.vars = c("x_index"))
sample.num <- length(unique(data.melt\$x_index))
attr.num <- length(unique(data.melt\$variable))

if(isbracket == T) {
    data.melt\$value <- gsub("\\\\(.*?\\\\)", "", data.melt\$value)
    data.melt\$value <- gsub("\\\\s+", "", data.melt\$value)
    data.melt\$value <- as.numeric(data.melt\$value)
}

barwidth = 0.8
y_expand = c(0.03, 0.03)
if(sample.num <= 10) {
    barwidth = 0.6
    width = 8
} else if(sample.num <= 30) {
    width = 11
    y_expand = c(0.02, 0.02)
} else if(sample.num <= 100) {
    width = 11 + 8 * (sample.num - 30) / 70
    y_expand = c(0.015, 0.015)
} else if(sample.num <= 200) {
    width = 20
    y_expand = c(0.012, 0.012)
} else {
    width = 30
    y_expand = c(0.002, 0.002)
}
if(attr.num > 1 && islegend == T) {
    legend.name <- as.character(unique(data.melt\$variable))
    legend.name.len <- max(nchar(legend.name))
    legend.width <- 1 * ceiling(attr.num / 20) * (legend.name.len / 20)
    width <- width + legend.width 
}
if(width > 20) {
    width <- 20
}

axis.size = AxisTextSize(data.melt\$x_index, width)
text.size = GeomTextSize(data.melt\$value, sample.num = sample.num, width = width)

y_min <- min(data.melt\$value)
y_max <- max(data.melt\$value) * 1.12
if(y_min > 0) {
    y_min = 0
}

RSC

    if($self->{pic}->{mulcolor} =~ /^T$|^True$/i and $attr_cnt == 1) {
        $rscript .= <<RSC;
if(length(colors) > sample.num) {
    colors <- colors[1:sample.num]
} else if(length(colors) != sample.num) {
    colors <- colorRampPalette(c(colors))(sample.num)
}
p <- ggplot(data.melt, aes(x = x_index, y = value, fill = x_index, color = I("#FFFFFF")))
RSC
    } else {
        $rscript .= <<RSC;
if(attr.num == 1) islegend = F
if(length(colors) > attr.num) {
    colors <- colors[1:attr.num]
} else if(length(colors) != attr.num) {
    colors <- colorRampPalette(c(colors))(attr.num)
}
p <- ggplot(data.melt, aes(x = x_index, y = value, fill = variable, color = I("#FFFFFF")))
RSC
    }

    if($self->{pic}->{bardodge} =~ /^T$|^True$/i and $attr_cnt > 1) {
        $rscript .= <<RSC;
p <- p + geom_bar(stat = "identity", position = position_dodge(), width = barwidth, size = 0.1)
RSC
    } elsif ($self->{pic}->{barfill} =~ /^T$|^True$/i and $attr_cnt > 1) {
        $rscript .= <<RSC;
p <- p + geom_bar(stat = "identity", position = position_fill(), width = barwidth, size = 0.1)
RSC
    } else {
        $rscript .= <<RSC;
p <- p + geom_bar(stat = "identity", position = position_stack(), width = barwidth, size = 0.1)
RSC
    }

    if($self->{pic}->{label} =~ /^T$|^True$/i and $attr_cnt == 1) {
        $rscript .= <<RSC;
p <- p + geom_text(aes(label = value), color = "#000000", size = text.size, 
    hjust = 0.5, vjust = -0.5) + scale_y_continuous(expand = y_expand, limits = c(y_min, y_max), 
    breaks = pretty_breaks(5))
RSC
    } elsif ($self->{pic}->{label} =~ /^T$|^True$/i and $self->{pic}->{bardodge} =~ /^T$|^True$/i
        and $attr_cnt > 1) {
        $rscript .= <<RSC;
text.size = GeomTextSize(data.melt\$value, sample.num = sample.num * attr.num, width = width)
p <- p + geom_text(aes(label = value), color = "#000000", size = text.size, 
    position = position_dodge(width = barwidth), vjust = -0.5, hjust = 0.5) + 
scale_y_continuous(expand = y_expand, breaks = pretty_breaks(5), limits = c(y_min, y_max))
RSC
    } else {
        $rscript .= <<RSC;
p <- p + scale_y_continuous(expand = y_expand, breaks = pretty_breaks(5))
RSC
    }

    $rscript .= <<RSC;

p <- p + scale_fill_manual(values = colors) +
labs(x = xlab, y = ylab, title = title) +
mytheme +
theme(axis.text.x = element_text(size = axis.size)) + 
xlab_angle(unique(data.melt\$x_index), width = width)

if(islegend != T) {
    width <- width - 1
    p <- p + theme(legend.position = "none")
}

if(iscoordflip == T) {
    p <- p + coord_flip() + 
    theme(axis.text.y = element_text(size = axis.size), 
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5, vjust = 0.5))
    tmp_width <- width
    tmp_height <- height
    width <- tmp_height
    height <- tmp_width
}

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "bar";
    $self->{obj} = 1;

    return $self;

}

sub error_bar {

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

    $self->{pic}->{header}   = $opts{'-header'}   ? $opts{'-header'}   : "T";
    $self->{pic}->{mulcolor} = $opts{'-mulcolor'} ? $opts{'-mulcolor'} : "F";
    $self->{pic}->{format}   = $opts{'-format'}   ? $opts{'-format'}   : "matrix";  # matrix | list3

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";    # 'red,blue'

    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{label}   = $opts{'-label'}  ? $opts{'-label'}  : "F";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 6;

    # color
    my $color = "";
    my $label_cnt  = 0;
    my $sample_cnt = 0;
    my @return = data_check($self->{pic}->{file}, $self->{pic}->{format});
    if($self->{pic}->{format} eq 'matrix') {
        #| for matrix format, we assume that:
        #| sample_id    bar_height    error_bar_height
        #| ...
        #| or :
        #| sample_id    bar_height    error_bar_height    label
        #| ...
        $label_cnt  = $return[0];
        $sample_cnt = $return[0];
    }
    elsif($self->{pic}->{format} eq 'list3') {
        #| not finished yet
        #| for list3 format, we assume that:
        #| id1   id2    id3
        #| v1    v2     v3
        #| v4
        $label_cnt  = $return[0];
        $sample_cnt = $return[0];
    }
    else {
        ERROR("error_bar", "the -format [$self->{pic}->{format}] is wrong. ");
    }
    if($self->{pic}->{header} =~ /^T$|^True$/i) {
        $label_cnt --;
        $sample_cnt --;
    }
    if($self->{pic}->{mulcolor} !~ /^T$|^True$/i) {
        $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "F";
        $label_cnt = 1;
    }

    if($label_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("bar", "default", $label_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $label_cnt) {
            $color = join(",", map { "'$_'" } @$colors[0 .. ($label_cnt - 1)]);
        }
        else {
            $color = join(",", map { "'$_'" } @$colors);
            $color = "colorRampPalette(c($color))($label_cnt)";
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
    elsif($self->{pic}->{collist} ne "none") {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            $color = join(",", map {"'$_'"} @colors);
        }
        else {
            $color = "'$self->{pic}->{collist}'";
        }
    }

    # theme
    my $header = ggplot2_head();
    my $theme  = ggplot2_bar_theme_default();

    my $rscript = <<RSC;
$header
$theme
source("$SOFT{'ggplot2'}")

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

title   <- "$self->{pic}->{title}"
xlab    <- "$self->{pic}->{xlab}"
ylab    <- "$self->{pic}->{ylab}"
islabel   <- as.logical("$self->{pic}->{label}")
islegend  <- as.logical("$self->{pic}->{legend}")

width   <- $self->{pic}->{width}
height  <- $self->{pic}->{height}
colors <- c($color)

RSC

    if($self->{pic}->{format} eq 'matrix') {
        $rscript .= <<RSC;
data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
if(length(data) == 4) {
    colnames(data) = c("x_index", "bar_height", "errorbar_height", "label")
} else if(length(data) == 3) {
    colnames(data) = c("x_index", "bar_height", "errorbar_height")
}
sample.num <- length(unique(data\$x_index))
if(length(colors) > sample.num) {
    colors <- colors[1:sample.num]
} else if(length(colors) != sample.num) {
    colors <- colorRampPalette(c(colors))(sample.num)
}

barwidth = 0.8
y_expand = c(0.03, 0.03)
if(sample.num <= 10) {
    barwidth = 0.6
    width = 7
} else if(sample.num <= 30) {
    width = 10
    y_expand = c(0.02, 0.02)
} else if(sample.num <= 100) {
    width = 10 + 8 * (sample.num - 30) / 70
    y_expand = c(0.02, 0.02)
} else {
    width = 20
    y_expand = c(0.015, 0.015)
}

axis.size = AxisTextSize(data\$x_index, width)
text.size = 4
if(ncol(data) == 4) {
    text.size = GeomTextSize(data\$x_index, sample.num = sample.num, width = width)
}

y_min <- min(data\$bar_height) 
y_max <- max(data\$bar_height)
ye_min <- min(data\$bar_height - data\$errorbar_height) * 1.12
ye_max <- max(data\$bar_height + data\$errorbar_height) * 1.12
if(ye_min > 0) {
    ye_min = 0
}

p <- ggplot(data, aes(x = x_index, y = bar_height, fill = x_index, color = I("#FFFFFF"))) 
if(y_max > 0) {
    p <- p + 
    geom_errorbar(aes(ymin = bar_height, ymax = bar_height + errorbar_height), 
        color = "#000000", position = position_dodge(width = barwidth), 
        width = barwidth * 0.6, size = 0.6) +
    scale_y_continuous(expand = y_expand, limits = c(ye_min, ye_max))
} else {
    p <- p + 
    geom_errorbar(aes(ymin = bar_height, ymax = bar_height - errorbar_height), 
        color = "#000000", position = position_dodge(width = barwidth), 
        width = barwidth * 0.6, size = 0.6) +
    scale_y_continuous(expand = y_expand, limits = c(ye_min, ye_max))
}
p <- p + geom_bar(stat = "identity", position = position_dodge(), width = barwidth, size = 0.6)

if(length(data) == 4 && islabel == T) {
    if(y_min > 0) {
        p <- p + geom_text(aes(y = y_min * 0.2, label = label), color = "#000000", 
        size = text.size, hjust = 0.5, vjust = 0.5)
    } else {
        p <- p + geom_text(aes(y = y_min * 0.1, label = label), color = "#000000", 
        size = text.size, hjust = 0.5, vjust = 0.5)
    }
}

if(y_min < 0) {
    p <- p + geom_hline(yintercept = 0)
}

p <- p + scale_fill_manual(values = colors) +
labs(x = xlab, y = ylab, title = title) +
mytheme +
theme(axis.text.x = element_text(size = axis.size)) +
xlab_angle(unique(data\$x_index), width = width)

if(islegend == F) {
    p <- p + theme(legend.position = "none")
}

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC
    }
    elsif($self->{pic}->{format} eq 'list') {
        $rscript = <<RSC;

RSC
    }

    $self->{rscript} = $rscript;
    $self->{method} = "errorbar";
    $self->{obj} = 1;

    return $self;
}

sub polar_bar {

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

    $self->{pic}->{collist}     = defined $opts{'-collist'}     ? $opts{'-collist'} : "none";
    $self->{pic}->{colset}      = defined $opts{'-colset'}      ? $opts{'-colset'} : "rainbow";
    $self->{pic}->{title}       = defined $opts{'-title'}       ? $opts{'-title'} : "none";
    $self->{pic}->{format}      = defined $opts{'-format'}      ? $opts{'-format'} : "none";
    $self->{pic}->{header}      = defined $opts{'-header'}      ? $opts{'-header'} : "T";
    $self->{pic}->{show_border} = defined $opts{'-show_border'} ? $opts{'-show_border'} : "F";
    $self->{pic}->{show_label}  = defined $opts{'-show_label'}  ? $opts{'-show_label'}  : "T";
    $self->{pic}->{show_count}  = defined $opts{'-show_count'}  ? $opts{'-show_count'}  : "T";
    $self->{pic}->{direction}   = defined $opts{'-direction'}   ? $opts{'-direction'} : "clockwise";
    $self->{pic}->{display}     = defined $opts{'-display'}     ? $opts{'-display'} : "type1";
    $self->{pic}->{width}       = defined $opts{'-width'}       ? $opts{'-width'}  : 8;
    $self->{pic}->{height}      = defined $opts{'-height'}      ? $opts{'-height'} : 8;

    $self->{convert_density} = defined $opts{'-convert_density'} ? $opts{'-convert_density'} : 1000;
    $self->{pic}->{collist} = "'$self->{pic}->{collist}'";
    $self->{pic}->{title} = "'$self->{pic}->{title}'";

    data_check($self->{pic}->{file}, "polarbar");

    my $shell = "export LANG='zh_CN.UTF-8' ; mkdir -p $self->{pic}->{outdir} ; ". 
    "$SOFT{Rscript} $SOFT{polar_bar} $self->{pic}->{file} ".
    "$self->{pic}->{outdir}/$self->{pic}->{outname} ".
    "$self->{pic}->{collist} $self->{pic}->{colset} $self->{pic}->{title} $self->{pic}->{format} ".
    "$self->{pic}->{header} $self->{pic}->{show_border} $self->{pic}->{show_label} ".
    "$self->{pic}->{show_count} $self->{pic}->{direction} $self->{pic}->{display} ".
    "$self->{pic}->{width} $self->{pic}->{height}";

    $self->{rscript} = "";
    $self->{command} = $shell;
    $self->{method} = "polarbar";
    $self->{obj} = 1;

    return $self;

}

sub sankey {

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

    $self->{pic}->{collist}      = defined $opts{'-collist'}      ? $opts{'-collist'} : "none";
    $self->{pic}->{colset}       = defined $opts{'-colset'}       ? $opts{'-colset'} : "set1";
    $self->{pic}->{transparency} = defined $opts{'-transparency'} ? $opts{'-transparency'} : 0.7;
    $self->{pic}->{band_width}   = defined $opts{'-band_width'}   ? $opts{'-band_width'} : 0.15;
    $self->{pic}->{band_knotpos} = defined $opts{'-band_knotpos'} ? $opts{'-band_knotpos'} : 0.4;
    $self->{pic}->{title}        = defined $opts{'-title'}        ? $opts{'-title'} : "none";
    $self->{pic}->{xlab}         = defined $opts{'-xlab'}         ? $opts{'-xlab'} : "none";
    $self->{pic}->{ylab}         = defined $opts{'-ylab'}         ? $opts{'-ylab'} : "none";
    $self->{pic}->{title_size}   = defined $opts{'-title_size'}   ? $opts{'-title_size'} : 16;
    $self->{pic}->{axis_size}    = defined $opts{'-axis_size'}    ? $opts{'-axis_size'} : 12;
    $self->{pic}->{label_size}   = defined $opts{'-label_size'}   ? $opts{'-label_size'} : 4;
    $self->{pic}->{label_color}  = defined $opts{'-label_color'}  ? $opts{'-label_color'} : "black";
    $self->{pic}->{fill_column}  = defined $opts{'-fill_column'}  ? $opts{'-fill_column'} : "default";
    $self->{pic}->{stratum_color} = defined $opts{'-stratum_color'} ? $opts{'-stratum_color'} : "black";
    $self->{pic}->{alluvium_color} = defined $opts{'-alluvium_color'} ? $opts{'-alluvium_color'} : "NA";

    $self->{pic}->{collist} = "'$self->{pic}->{collist}'";
    $self->{pic}->{colset} = "'$self->{pic}->{colset}'";
    $self->{pic}->{title} = "'$self->{pic}->{title}'";
    $self->{pic}->{xlab} = "'$self->{pic}->{xlab}'";
    $self->{pic}->{ylab} = "'$self->{pic}->{ylab}'";

    data_check($self->{pic}->{file}, 'sankey');

    my $shell = "export LANG='zh_CN.UTF-8' ; mkdir -p $self->{pic}->{outdir} ; ". 
    "$SOFT{Rscript} $SOFT{sankey_plot} $self->{pic}->{file} ".
    "$self->{pic}->{outdir}/$self->{pic}->{outname} ".
    "$self->{pic}->{collist} $self->{pic}->{colset} $self->{pic}->{transparency} ".
    "$self->{pic}->{band_width} $self->{pic}->{band_knotpos} ".
    "$self->{pic}->{title} $self->{pic}->{xlab} $self->{pic}->{ylab} ".
    "$self->{pic}->{title_size} $self->{pic}->{axis_size} $self->{pic}->{label_size} ".
    "$self->{pic}->{label_color} $self->{pic}->{stratum_color} ".
    "$self->{pic}->{alluvium_color} $self->{pic}->{fill_column} ";

    $self->{rscript} = "";
    $self->{command} = $shell;
    $self->{method} = "sankey";
    $self->{obj} = 1;

    return $self;

}

return 1;

