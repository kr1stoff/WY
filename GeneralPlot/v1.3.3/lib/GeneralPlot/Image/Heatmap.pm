package GeneralPlot::Image::Heatmap;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.10 16:13:50          |
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
        'type'   => 'heatmap',
        'method' => 'none',
        'obj'    => 0,
        'save'   => 0,
        'run'    => 0,
    };

    bless $self, $class;

    return $self;

}

#-------------------------------------------------------------------------------
#  general a heatmap 
#-------------------------------------------------------------------------------
sub plot {

    my ($self, %opts) = @_;

    $opts{'-method'} ||= "none";

    if(defined $opts{'-method'} and $opts{'-method'} eq 'corr') {
        $self->{method} = "corr";
        heatmap_corr($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'ggcor') {
        $self->{method} = "ggcor";
        heatmap_ggcor($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'cluster') {
        $self->{method} = "cluster";
        heatmap_cluster($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'heatmap') {
        $self->{method} = "heatmap";
        heatmap_base($self, %opts);
    }
    else {
        ERROR("$self->{type}", "check the -method [$opts{'-method'}]. ");
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general a default setting heatmap (not finished) 
#-------------------------------------------------------------------------------
sub heatmap_base {

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

    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{rowcol}  = $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";     # row
    $self->{pic}->{namefix} = $opts{'-namefix'} ? $opts{'-namefix'} : "none";    # _fpkm ...
    $self->{pic}->{corr}    = $opts{'-corr'}    ? $opts{'-corr'}    : "pearson"; # kendall | spearman

    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";  # heatmap_default_4
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";  # 'red,blue'

    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{label}   = $opts{'-label'}  ? $opts{'-label'}  : "F";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 7;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 7;

    # color settings
    my $color = "";
    my $colors = fetch_color(-type => "heatmap", -num => 0);
    $color = join(",", map {"'$_'"} @$colors);
    if($self->{pic}->{colset} ne 'none') {
        my @set = split /\_/, $self->{pic}->{colset};
        if(@set == 2) {
            $colors = fetch_color(-type => $set[0], -tag => $set[1], -num => 0);
        }
        elsif(@set == 3) {
            $colors = fetch_color(-type => $set[0], -tag => $set[1], -num => $set[2]);
        }
        else {
            WARN("unknown color set. ");
        }
        $color = join(",", map {"'$_'"} @$colors);
    }

    if($self->{pic}->{collist} ne 'none') {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            $color = join(",", map {"'$_'"} @colors);
        }
    }

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_heatmap_theme_default();

    my $rscript = <<RSC;
$header
$theme

source("$SOFT{'ggplot2'}")

readTable <- function (file, rowcol) {
    data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, 
    row.names = 1, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

RSC

    if(!defined $self->{stat}) {
        $rscript .= <<RSC;
file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

RSC
    }
    else {
        $rscript .= <<RSC;
$self->{stat}

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

RSC
    }

    $rscript .= <<RSC;
title  <- "$self->{pic}->{title}"
xlab   <- "$self->{pic}->{xlab}"
ylab   <- "$self->{pic}->{ylab}"
islabel  <- as.logical("$self->{pic}->{label}")
islegend <- as.logical("$self->{pic}->{legend}")

width  <- $self->{pic}->{width}
height <- $self->{pic}->{height}
colors <- c($color)

data <- readTable(file, "$self->{pic}->{rowcol}")
RSC

    if($self->{pic}->{namefix} ne 'none') {
        $rscript .= <<RSC;
getcol <- grep("$self->{pic}->{namefix}", colnames(data))
if(length(getcol) > 0) {
    data <- data[, getcol]
    colnames(data) <- gsub("$self->{pic}->{namefix}", "", colnames(data))
}
RSC
    }


    $rscript .= <<RSC;
data.cor <- cor(data, use = "p", method = "$self->{pic}->{corr}")
write.table(data.frame("Sample" = rownames(data.cor), data.cor, check.names = F), 
    file = paste0(outdir, "/", outname, ".correlation.xls"), sep = "\\t", quote = F, row.names = F)

data.cor.m <- melt(data.cor)
data.cor.m\$Var2 <- factor(data.cor.m\$Var2, levels = rev(unique(data.cor.m\$Var2)))
data.cor.m\$format_value <- sprintf(data.cor.m\$value, fmt = "%0.3f")
sample.num <- length(unique(data.cor.m\$Var2))

if(sample.num <= 10) {
    width = 8
    height = 8
} else if(sample.num <= 30) {
    width = 10
    height = 10
} else if(sample.num <= 100) {
    width = 10 + 8 * (sample.num - 30) / 70
    height = width
} else {
    width = 20
    height = 20
}

axis.size <- AxisTextSize(data.cor.m\$Var2, width = width)
text.size <- GeomTextSize(data.cor.m\$format_value, sample.num = sample.num, 
    width = width, legend = F)

p <- ggplot(data.cor.m, aes(Var1, Var2, fill = value)) + 
geom_tile(color = "white") +
scale_fill_gradientn(colors = colors) +
scale_y_discrete(position = "right") + 
labs(x = xlab, y = ylab, title = title) + 
coord_fixed(ratio = 1) +
mytheme

if(islabel == T) {
    p <- p + geom_text(aes(label = format_value), size = text.size, color = "grey10")
}
if(islegend != T) {
    p <- p + theme(legend.position = 'none')
}

p <- p + theme(axis.text = element_text(size = axis.size)) + 
guides(fill = guide_colorbar(barwidth = 0.8, barheight = 10))

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "corr";
    $self->{obj} = 1;

    return $self;

}



#-------------------------------------------------------------------------------
#  general the correlation heatmap 
#-------------------------------------------------------------------------------
sub heatmap_corr {

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

    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{rowcol}  = $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";     # row
    $self->{pic}->{namefix} = $opts{'-namefix'} ? $opts{'-namefix'} : "none";    # _fpkm ...
    $self->{pic}->{corr}    = $opts{'-corr'}    ? $opts{'-corr'}    : "pearson"; # kendall | spearman

    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";   # 'red,blue'

    $self->{pic}->{title}   = defined $opts{'-title'} ? $opts{'-title'} : "Sample Correlation";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{label}   = $opts{'-label'}  ? $opts{'-label'}  : "F";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 7;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 7;

    # color settings
    my $color = "";
    my ($type, $tag, $num) = ("heatmap", "corr", 0);
    if($self->{pic}->{colset} ne 'none') {
        my @set = split /\_/, $self->{pic}->{colset};
        $type = $set[0];
        $tag  = $set[1];
        $num  = defined $set[2] ? $set[2] : $num;
    }

    my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
    $color = join(",", map {"'$_'"} @$colors);

    if($self->{pic}->{collist} ne 'none') {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            $color = join(",", map {"'$_'"} @colors);
        }
    }

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_heatmap_theme_default();

    my $rscript = <<RSC;
$header
$theme

source("$SOFT{'ggplot2'}")

readTable <- function (file, rowcol) {
    data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, 
    row.names = 1, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

title  <- "$self->{pic}->{title}"
xlab   <- "$self->{pic}->{xlab}"
ylab   <- "$self->{pic}->{ylab}"
islabel  <- as.logical("$self->{pic}->{label}")
islegend <- as.logical("$self->{pic}->{legend}")

width  <- $self->{pic}->{width}
height <- $self->{pic}->{height}
colors <- c($color)

data <- readTable(file, "$self->{pic}->{rowcol}")
RSC

    if($self->{pic}->{namefix} ne 'none') {
        $rscript .= <<RSC;
getcol <- grep("$self->{pic}->{namefix}", colnames(data))
if(length(getcol) > 0) {
    data <- data[, getcol]
    colnames(data) <- gsub("$self->{pic}->{namefix}", "", colnames(data))
}
RSC
    }


    $rscript .= <<RSC;
data.cor <- cor(data, use = "p", method = "$self->{pic}->{corr}")
write.table(data.frame("Sample" = rownames(data.cor), data.cor, check.names = F), 
    file = paste0(outdir, "/", outname, ".correlation.xls"), sep = "\\t", quote = F, row.names = F)

data.cor.m <- melt(data.cor)
data.cor.m\$Var2 <- factor(data.cor.m\$Var2, levels = rev(unique(data.cor.m\$Var2)))
data.cor.m\$format_value <- sprintf(data.cor.m\$value, fmt = "%0.3f")
sample.num <- length(unique(data.cor.m\$Var2))

if(sample.num <= 10) {
    width = 8
    height = 8
} else if(sample.num <= 30) {
    width = 10
    height = 10
} else if(sample.num <= 100) {
    width = 10 + 8 * (sample.num - 30) / 70
    height = width
} else {
    width = 20
    height = 20
}

axis.size <- AxisTextSize(data.cor.m\$Var2, width = width)
text.size <- GeomTextSize(data.cor.m\$format_value, sample.num = sample.num, 
    width = width, legend = F)

p <- ggplot(data.cor.m, aes(Var1, Var2, fill = value)) + 
geom_tile(color = "white") +
scale_fill_gradientn(colors = colors) +
scale_y_discrete(position = "right") + 
labs(x = xlab, y = ylab, title = title) + 
coord_fixed(ratio = 1) +
mytheme

if(islabel == T) {
    p <- p + geom_text(aes(label = format_value), size = text.size, color = "grey10")
}
if(islegend != T) {
    p <- p + theme(legend.position = 'none')
}

p <- p + theme(axis.text = element_text(size = axis.size)) + 
guides(fill = guide_colorbar(barwidth = 0.8, barheight = 10))

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "corr";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  general the correlation heatmap with ggcor package
#-------------------------------------------------------------------------------
sub heatmap_ggcor {

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

    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{rowcol}  = $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";     # row
    $self->{pic}->{namefix} = $opts{'-namefix'} ? $opts{'-namefix'} : "none";    # _fpkm ...

    $self->{pic}->{cor}     = $opts{'-cor'}     ? $opts{'-cor'}     : "pearson"; # kendall | spearman
    $self->{pic}->{cortest} = $opts{'-cortest'} ? $opts{'-cortest'} : "F";
    $self->{pic}->{pvalue}  = $opts{'-pvalue'}  ? $opts{'-pvalue'}  : 0.05;
    $self->{pic}->{postype} = $opts{'-postype'} ? $opts{'-postype'} : "full";    # upper | lower

    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";   # 'red,blue'

    $self->{pic}->{title}   = defined $opts{'-title'} ? $opts{'-title'} : "Sample Correlation";
    $self->{pic}->{xlab}    = $opts{'-xlab'}    ? $opts{'-xlab'}    : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}    ? $opts{'-ylab'}    : "";
    $self->{pic}->{label}   = $opts{'-label'}   ? $opts{'-label'}   : "F";
    $self->{pic}->{legend}  = $opts{'-legend'}  ? $opts{'-legend'}  : "T";

    $self->{pic}->{width}   = $opts{'-width'}   ? $opts{'-width'}   : 7;
    $self->{pic}->{height}  = $opts{'-height'}  ? $opts{'-height'}  : 7;

    # color settings
    my $color = "";
    my ($type, $tag, $num) = ("heatmap", "corr", 0);
    if($self->{pic}->{colset} ne 'none') {
        my @set = split /\_/, $self->{pic}->{colset};
        $type = $set[0];
        $tag  = $set[1];
        $num  = defined $set[2] ? $set[2] : $num;
    }

    my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
    $color = join(",", map {"'$_'"} @$colors);

    if($self->{pic}->{collist} ne 'none') {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            $color = join(",", map {"'$_'"} @colors);
        }
    }

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_heatmap_theme_ggcor();

    my $rscript = <<RSC;
$header
library(ggcor)
$theme

source("$SOFT{'ggplot2'}")

readTable <- function (file, rowcol) {
    data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, 
    row.names = 1, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

RSC

    if(!defined $self->{stat}) {
        $rscript .= <<RSC;
file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

RSC
    }
    else {
        $rscript .= <<RSC;
$self->{stat}

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

RSC
    }

    $rscript .= <<RSC;
title  <- "$self->{pic}->{title}"
xlab   <- "$self->{pic}->{xlab}"
ylab   <- "$self->{pic}->{ylab}"
islabel  <- as.logical("$self->{pic}->{label}")
islegend <- as.logical("$self->{pic}->{legend}")

width  <- $self->{pic}->{width}
height <- $self->{pic}->{height}
colors <- c($color)

data <- readTable(file, "$self->{pic}->{rowcol}")
data.cor <- correlate(data, method = "$self->{pic}->{cor}", cor.test = $self->{pic}->{cortest})
data.cor.tb <- fortify_cor(data, method = "$self->{pic}->{cor}", type = "$self->{pic}->{postype}", 
    cor.test = $self->{pic}->{cortest})

write.table(data.frame("Sample" = rownames(data.cor\$r), data.cor\$r, check.names = F), 
    file = paste0(outdir, "/", outname, ".correlation.xls"), sep = "\\t", quote = F, row.names = F)

RSC

    if($self->{pic}->{cortest} =~ /^T$|^True$/i) {
        $rscript .= <<RSC;
p.value <- data.cor\$p.value
rownames(p.value) <- rownames(data.cor\$r)
colnames(p.value) <- colnames(data.cor\$r)
write.table(data.frame("Sample" = rownames(p.value), p.value, check.names = F), 
    file = paste0(outdir, "/", outname, ".pvalue.xls"), sep = "\\t", quote = F, row.names = F)

# filter data here, geom_cross() or geom_mark() may be better
data.cor.tb <- filter(data.cor.tb, p.value < $self->{pic}->{pvalue})

RSC
    }

    $rscript .= <<RSC;
sample.num <- length(unique(data.cor.tb\$".col.names"))

if(sample.num <= 10) {
    width = 8
    height = 8
} else if(sample.num <= 30) {
    width = 10
    height = 10
} else if(sample.num <= 100) {
    width = 10 + 8 * (sample.num - 30) / 70
    height = width
} else {
    width = 20
    height = 20
}

axis.size <- AxisTextSize(data.cor.tb\$".col.names", width = width)

p <- ggcor(data.cor.tb) +
geom_square(aes(r0 = r, fill = r)) +
scale_fill_gradientn(colors = colors) +
coord_fixed(expand = FALSE) +
geom_panel_grid(size = 1, color = "grey50") + 
geom_panel_grid(size = 0.6, color = "#FFFFFF") + 
labs(x = xlab, y = ylab, title = title) + 
mytheme + 
theme(axis.text = element_text(size = axis.size))

if(islabel == T) {
    cor.r.format <- sprintf("%0.3f", data.cor.tb\$r)
    text.size <- GeomTextSize(cor.r.format, sample.num = sample.num, width = width, legend = F)
    p <- p + geom_number(aes(num = r), size = text.size, color = "grey10", digits = 3)
}
if(islegend == F) {
    p <- p + theme(legend.position = 'none')
}

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "ggcor";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  general the cluster heatmap by pheatmap
#-------------------------------------------------------------------------------
sub heatmap_cluster {

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
    $self->{pic}->{rowcol}   = $opts{'-rowcol'}   ? $opts{'-rowcol'}   : "col";    # row
    $self->{pic}->{namefix}  = $opts{'-namefix'}  ? $opts{'-namefix'}  : "none";   # _fpkm ...
    $self->{pic}->{clustrow} = $opts{'-clustrow'} ? $opts{'-clustrow'} : "T";   # T | F
    $self->{pic}->{clustcol} = $opts{'-clustcol'} ? $opts{'-clustcol'} : "T";   # T | F

    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";   # 'red,blue'

    $self->{pic}->{title}   = $opts{'-title'}   ? $opts{'-title'}   : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}    ? $opts{'-xlab'}    : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}    ? $opts{'-ylab'}    : "";
    $self->{pic}->{label}   = $opts{'-label'}   ? $opts{'-label'}   : "F";
    $self->{pic}->{legend}  = $opts{'-legend'}  ? $opts{'-legend'}  : "T";

    $self->{pic}->{width}   = $opts{'-width'}   ? $opts{'-width'}   : 7;
    $self->{pic}->{height}  = $opts{'-height'}  ? $opts{'-height'}  : 7;

    # color settings
    my $color = "";
    my ($type, $tag, $num) = ("heatmap", "default", 0);
    if($self->{pic}->{colset} ne 'none') {
        my @set = split /\_/, $self->{pic}->{colset};
        $type = $set[0];
        $tag  = $set[1];
        $num  = defined $set[2] ? $set[2] : $num;
    }

    my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
    $color = join(",", map {"'$_'"} @$colors);

    if($self->{pic}->{collist} ne 'none' and $self->{pic}->{collist} =~ /,/) {
        my @colors = split /,/, $self->{pic}->{collist};
        $color = join(",", map {"'$_'"} @colors);
    }

    # theme settings
    # ...

    my $rscript = <<RSC;
library(pheatmap)

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

isscale        = 'row' ## scale or not
iscluster_row  = $self->{pic}->{clustrow}  ## cluster_row or not
iscluster_col  = $self->{pic}->{clustcol}  ## cluster_col or not
incell_width   = 'NA'  ## cell width 
incell_height  = 'NA'  ## cell height
isdisplay_num  = $self->{pic}->{label}     ## display number in cell
inborder_color = NA    ## border color
infontsize     = 10    ## fontsize
isshow_rowname = T     ## show rowname or not
isshow_colname = T     ## show rowname or not

readTable <- function (file, rowcol) {
    data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, 
    row.names = 1, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

## read data
data <- readTable(file, "$self->{pic}->{rowcol}")
RSC

    if($self->{pic}->{namefix} ne 'none') {
        $rscript .= <<RSC;
getcol <- grep("$self->{pic}->{namefix}", colnames(data))
if(length(getcol) > 0) {
    data <- data[, getcol]
    colnames(data) <- gsub("$self->{pic}->{namefix}", "", colnames(data))
}

RSC
    }

    $rscript .= <<RSC;
if(ncol(data) < 3) {
    isscale <- 'none'
    data[data == 0] <- 0.001
    data <- log10(data)
}
if(nrow(data) > 50) {
    isshow_rowname <- F
}

# delete the row while rowsum is 0
data <- data[which(rowSums(data) > 0),]

## deal color
colors <- colorRampPalette(c($color))(100)

## check nrow and ncol of data, 1 can't scale and cluster
if(nrow(data) == 1) {
    if(isscale == 'column') {
        cat("=>Attention: Can't scale by column with only one column!!", 
        "So will not scale the data!! \\n")
        isscale = 'none'
    }
    iscluster_row = 'FALSE'
}

if(ncol(data) == 1) {
    if(isscale == 'row') {
        cat("=>Attention: Can't scale by row with only one row!! So will not scale the data!! \\n")
        isscale = 'none'
    }
    iscluster_col = 'FALSE'
}

## too many rows, can't defined cell width and height
if(nrow(data) > 1900) {
    incell_width = 'NA'
    incell_height = 'NA'
}

if(nrow(data) > 10000) {
    cat("[pheatmap warnings] genes more than 10000, use the first 10000 only \\n")
    data = data[1:10000, ]
}

## get maxchar of rowname, if maxchar too big, the fig will be small
maxchar <- max(nchar(rownames(data), allowNA = TRUE), na.rm = TRUE)

if(maxchar > 40) {
    if(incell_width == 'NA') {
        incell_width = 25
    }
}
if(incell_width == 'NA') {
    incell_width = as.logical(incell_width)
}
if(incell_height == 'NA') {
    incell_height = as.logical(incell_height)
}

setwd(outdir)

outfile = paste0(outdir, "/", outname, ".pdf")
pheatmap(data,
    scale = isscale,  # normalization
    cluster_rows = as.logical(iscluster_row),  # cluster
    cluster_cols = as.logical(iscluster_col),  # cluster
    color = colors,  # color
    cellwidth = as.numeric(incell_width), 
    cellheight = as.numeric(incell_height),  #cell width and height
    display_numbers = isdisplay_num,  # display number in the cells
    border = as.character(inborder_color),  # cell border
    fontsize = infontsize,
    show_rownames = isshow_rowname,
    show_colnames = isshow_colname,
    filename = outfile  # save file
)

outfile = paste0(outdir, "/", outname, ".png")
p.msg <- pheatmap(data,
    scale = isscale,  # normalization
    cluster_rows = as.logical(iscluster_row),  # cluster
    cluster_cols = as.logical(iscluster_col),  # cluster
    color = colors,  # color
    cellwidth = as.numeric(incell_width), 
    cellheight = as.numeric(incell_height),  #cell width and height
    display_numbers = isdisplay_num,  # display number in the cells
    border = as.character(inborder_color),  # cell border
    fontsize = infontsize,
    show_rownames = isshow_rowname,
    show_colnames = isshow_colname,
    filename = outfile  # save file
)

files = dir()
if("Rplots.pdf" %in% files)  res = file.remove("Rplots.pdf")
q()

if(iscluster_row == "TRUE" || iscluster_col == "TRUE") {
    outfile = paste0(outdir, "/", "pheatmap_reorder.xls")
    order_row = 1:nrow(data)
    order_col = 1:ncol(data)
    if(iscluster_row == "TRUE"){
        order_row = p.msg\$tree_row\$order
    }
    if(iscluster_col == "TRUE"){
        order_col = p.msg\$tree_col\$order
    }
    
    datat = data.frame(data[order_row, order_col])
    datat = data.frame(rownames(datat), datat, check.names = F)
    colnames(datat)[1] = "id"
    write.table(datat, file = outfile, row.names = FALSE, quote = FALSE, sep = '\\t')
}

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "cluster";
    $self->{obj} = 1;

    return $self;

}

return 1;

