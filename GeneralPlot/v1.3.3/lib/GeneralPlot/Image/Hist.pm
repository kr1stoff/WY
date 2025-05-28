package GeneralPlot::Image::Hist;
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
#  new a object 
#-------------------------------------------------------------------------------
sub new {

    my ($class, %opts) = @_;

    my $self = {
        'type'   => 'hist',
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

    $opts{'-method'} ||= "none";

    if(defined $opts{'-method'} and $opts{'-method'} eq 'hist') {
        $self->{method} = $opts{'-method'};
        hist($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'density') {
        $self->{method} = $opts{'-method'};
        density($self, %opts);
    }
    else {
        ERROR("$self->{type}", "check the -method [$opts{'-method'}]. ");
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general a single hist chart from a format table directly
#-------------------------------------------------------------------------------
sub hist {

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

    $self->{pic}->{format}   = $opts{'-format'}   ? $opts{'-format'}   : "table";  # table | matrix
    $self->{pic}->{header}   = $opts{'-header'}   ? $opts{'-header'}   : "T";
    $self->{pic}->{colname}  = $opts{'-colname'}  ? $opts{'-colname'}  : "NULL";
    $self->{pic}->{facet}    = $opts{'-facet'}    ? $opts{'-facet'}    : "F";
    $self->{pic}->{identity} = $opts{'-identity'} ? $opts{'-identity'} : "F";
    $self->{pic}->{density}  = $opts{'-density'}  ? $opts{'-density'}  : "F";

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";    # 'red,blue'

    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 8;
    $self->{pic}->{bins}    = $opts{'-bins'}   ? $opts{'-bins'}   : 30;

    # color
    my ($sample_cnt, $data_type) = (0, 0);
    if($self->{pic}->{format} eq 'table') {
        ($sample_cnt, $data_type) = count_table($self->{pic}->{file}, $self->{pic}->{header});
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix");
        $sample_cnt = $ncol - 1;
        $data_type = 2;
    }
    else {
        ERROR("hist", "check the -format [$self->{pic}->{format}]. ");
    }

    my $color = "";
    if($sample_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("hist", "default", $sample_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $sample_cnt) {
            $color = join(",", map { "'$_'" } @$colors[0 .. ($sample_cnt - 1)]);
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
    my $theme  = ggplot2_hist_theme_default();

    my $rscript = <<RSC;
$header
$theme

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"
colname <- "$self->{pic}->{colname}"

RSC

    if($self->{pic}->{format} eq 'table') {
        #| the table format is like:
        #| value
        #| ...
        #| or :
        #| sample    value
        #| ...

        $rscript .= <<RSC;
data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
if(ncol(data) == 1) {
    colnames(data)[1] <- "V_alue"
} else if(ncol(data) == 2) {
    colnames(data) <- c("S_ample", "V_alue")
} else {
    if(colname == "NULL") {
        colnames(data) <- c("S_ample", "V_alue")
    } else {
        colnames(data)[1] <- "S_ample"
        colnames(data)[colnames(data) == colname] <- "V_alue"
    }
}

RSC
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        #| the matrix format is like:
        #| attr    sample1    sample2    ...
        #| attr1   value1     value2     ...
        #| ...

        $rscript .= <<RSC;
data.raw <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
colnames(data.raw)[1] <- "X_index"
data.melt <- melt(data.raw, id.vars = c("X_index"))
data <- data.frame("S_ample" = data.melt\$variable, "V_alue" = data.melt\$value)

RSC
    }

    $rscript .= <<RSC;
title    <- "$self->{pic}->{title}"
xlab     <- "$self->{pic}->{xlab}"
ylab     <- "$self->{pic}->{ylab}"
islegend <- as.logical("$self->{pic}->{legend}")

width = $self->{pic}->{width}
height = $self->{pic}->{height}
colors <- c($color)

bins = $self->{pic}->{bins}
yexpand = c(0.02, 0.00)

y_min <- min(data\$V_alue)
y_max <- max(data\$V_alue)

if(length(data\$S_ample) == 0) {
    colors <- colors[1]
    p <- ggplot(data, aes(x = V_alue, fill = I(colors), color = I("#FFFFFF"))) 
} else {
    sample.num <- length(unique(data\$S_ample))    
    if(length(colors) > sample.num) {
        colors <- colors[1:sample.num]
    } else if(length(colors) != sample.num) {
        colors <- colorRampPalette(c(colors))(sample.num)
    }

    data\$S_ample <- factor(data\$S_ample, levels = unique(data\$S_ample), order = T)
    p <- ggplot(data, aes(x = V_alue, fill = S_ample, color = I("#FFFFFF"))) 
}

RSC

    if($self->{pic}->{density} =~ /^T$|^True$/i) {
        $rscript .= <<RSC;
p <- p + geom_histogram(aes(y = ..density..), position = position_stack(), 
    bins = bins, size = 0.1) + 
geom_density(color = "#000000", fill = NA, alpha = 0.5) + 
scale_x_continuous(breaks = pretty_breaks(5)) + 
scale_y_continuous(breaks = pretty_breaks(5), expand = yexpand)

RSC
    }
    elsif($self->{pic}->{identity} =~ /^T$|^True$/i and $data_type != 1) {
        $rscript .= <<RSC;
if(length(data\$S_ample) != 0) {
    p <- p + geom_histogram(position = position_identity(), bins = bins, 
        size = 0.1, alpha = 0.5) +
    scale_x_continuous(breaks = pretty_breaks(5)) + 
    scale_y_continuous(breaks = pretty_breaks(5), expand = yexpand)
}

RSC
    }
    else {
        $rscript .= <<RSC;
p <- p + geom_histogram(position = position_stack(), bins = bins, size = 0.1) + 
scale_x_continuous(breaks = pretty_breaks(5)) + 
scale_y_continuous(breaks = pretty_breaks(5), expand = yexpand)

RSC
    }

    if($self->{pic}->{facet} =~ /^T$|^True$/i and $data_type != 1) {
        $rscript .= <<RSC;
if(length(data\$S_ample) != 0) {
    p <- p + facet_wrap(.~S_ample)
}

RSC
    }
    elsif($self->{pic}->{density} =~ /^T$|^True$/i and $data_type != 1) {
        $rscript .= <<RSC;
if(length(data\$S_ample) != 0) {
    p <- p + facet_wrap(.~S_ample)
}

RSC
    }

    $rscript .= <<RSC;
if(length(data\$S_ample) != 0) {
    p <- p + scale_fill_manual(values = colors)
}
p <- p + labs(x = xlab, y = ylab, title = title) + mytheme 

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "hist";
    $self->{obj} = 1;

    return $self;

}

sub density {

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

    $self->{pic}->{format}   = $opts{'-format'}   ? $opts{'-format'}   : "table";  # table | matrix
    $self->{pic}->{header}   = $opts{'-header'}   ? $opts{'-header'}   : "T";
    $self->{pic}->{colname}  = $opts{'-colname'}  ? $opts{'-colname'}  : "NULL";
    $self->{pic}->{facet}    = $opts{'-facet'}    ? $opts{'-facet'}    : "F";
    $self->{pic}->{inline}   = $opts{'-inline'}   ? $opts{'-inline'}   : "F";

    $self->{pic}->{namefix} = defined $opts{'-namefix'} ? $opts{'-namefix'} : "none"; # _fpkm ...
    $self->{pic}->{filter0} = defined $opts{'-filter0'} ? $opts{'-filter0'} : "F";    # T | F
    $self->{pic}->{scale}   = defined $opts{'-scale'}   ? $opts{'-scale'}   : "none"; # log2 | log10 

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";    # 'red,blue'

    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 7;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 7;

    # color
    my ($sample_cnt, $data_type) = (0, 0);
    if($self->{pic}->{format} eq 'table') {
        ($sample_cnt, $data_type) = count_table($self->{pic}->{file}, $self->{pic}->{header});
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        if($self->{pic}->{namefix} eq 'none') {
            my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix");
            $sample_cnt = $ncol - 1;
            $data_type = 2;
        }
        else {
            my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix", $self->{pic}->{namefix});
            $sample_cnt = $ncol;
            $data_type = 2;
        }
    }
    else {
        ERROR("hist", "check the -format [$self->{pic}->{format}]. ");
    }

    my $color = "";
    if($sample_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("hist", "default", $sample_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $sample_cnt) {
            $color = join(",", map { "'$_'" } @$colors[0 .. ($sample_cnt - 1)]);
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
    my $theme  = ggplot2_hist_theme_default();

    my $rscript = <<RSC;
$header
$theme

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"
colname <- "$self->{pic}->{colname}"
scale   <- "$self->{pic}->{scale}"
filter0 <- as.logical($self->{pic}->{filter0})
namefix <- "$self->{pic}->{namefix}"

RSC

    if($self->{pic}->{format} eq 'table') {
        #| the table format is like:
        #| value
        #| ...
        #| or :
        #| sample    value
        #| ...

        $rscript .= <<RSC;
data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
if(ncol(data) == 1) {
    colnames(data)[1] <- "V_alue"
} else if(ncol(data) == 2) {
    colnames(data) <- c("S_ample", "V_alue")
} else {
    if(colname == "NULL") {
        colnames(data) <- c("S_ample", "V_alue")
    } else {
        colnames(data)[1] <- "S_ample"
        colnames(data)[colnames(data) == colname] <- "V_alue"
    }
}

RSC
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        #| the matrix format is like:
        #| attr    sample1    sample2    ...
        #| attr1   value1     value2     ...
        #| ...

        $rscript .= <<RSC;
data.raw <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
if(namefix != "none") {
    getcol <- grep(namefix, colnames(data.raw))
    if(length(getcol) > 0) {
        data.raw <- data.raw[, getcol]
        colnames(data.raw) <- gsub(namefix, "", colnames(data.raw))
    } else {
        stop("Error: no data found, check it. \\n")
    }
    colnames(data.raw)[1] <- "X_index"
    data.melt <- melt(data.raw, id.vars = c("X_index"))
    data <- data.frame("S_ample" = data.melt\$variable, "V_alue" = data.melt\$value)
} else {
    colnames(data.raw)[1] <- "X_index"
    data.melt <- melt(data.raw, id.vars = c("X_index"))
    data <- data.frame("S_ample" = data.melt\$variable, "V_alue" = data.melt\$value)
}

if(filter0 == T) {
    data <- data[data\$V_alue > 0, ]
}

# format the data
if(scale == "log2+1") {
    data\$V_alue <- log2(data\$V_alue + 1)
} else if(scale == "log10+1") {
    data\$V_alue <- log10(data\$V_alue + 1)
} else if(scale == "log2") {
    data\$V_alue[data\$V_alue == 0] = 0.001
    data\$V_alue <- log2(data\$V_alue)
} else if(scale == "log10") {
    data\$V_alue[data\$V_alue == 0] = 0.001
    data\$V_alue <- log10(data\$V_alue)
}

RSC
    }

    $rscript .= <<RSC;
title    <- "$self->{pic}->{title}"
xlab     <- "$self->{pic}->{xlab}"
ylab     <- "$self->{pic}->{ylab}"
islegend <- as.logical("$self->{pic}->{legend}")

width = $self->{pic}->{width}
height = $self->{pic}->{height}
colors <- c($color)

yexpand = c(0.02, 0.00)

RSC

    if($self->{pic}->{inline} =~ /^T$|^True$/i) {
        $rscript .= <<RSC;
if(length(data\$S_ample) == 0) {
    colors <- colors[1]
    p <- ggplot(data, aes(x = V_alue, color = I(colors))) 
} else {
    sample.num <- length(unique(data\$S_ample))    
    if(length(colors) > sample.num) {
        colors <- colors[1:sample.num]
    } else if(length(colors) != sample.num) {
        colors <- colorRampPalette(c(colors))(sample.num)
    }
    data\$S_ample <- factor(data\$S_ample, levels = unique(data\$S_ample), order = T)
    p <- ggplot(data, aes(x = V_alue, color = S_ample)) + 
    scale_color_manual(values = colors)
}

p <- p + geom_density(fill = NA) +
scale_x_continuous(breaks = pretty_breaks(5)) +
scale_y_continuous(breaks = pretty_breaks(5), expand = yexpand)

RSC

        if($self->{pic}->{facet} =~ /^T$|^True$/i and $data_type != 1) {
            $rscript .= <<RSC;
if(length(data\$S_ample) != 0) {
    p <- p + facet_wrap(.~S_ample)
}

RSC
        }
    }
    else {
        $rscript .= <<RSC;
if(length(data\$S_ample) == 0) {
    colors <- colors[1]
    p <- ggplot(data, aes(x = V_alue, fill = I(colors))) 
} else {
    sample.num <- length(unique(data\$S_ample))    
    if(length(colors) > sample.num) {
        colors <- colors[1:sample.num]
    } else if(length(colors) != sample.num) {
        colors <- colorRampPalette(c(colors))(sample.num)
    }
    data\$S_ample <- factor(data\$S_ample, levels = unique(data\$S_ample), order = T)
    p <- ggplot(data, aes(x = V_alue, fill = S_ample)) + 
    scale_fill_manual(values = colors)
}

RSC

        if($self->{pic}->{facet} =~ /^T$|^True$/i and $data_type != 1) {
            $rscript .= <<RSC;
p <- p + geom_density(color = "#FFFFFF", size = 0.2) +
scale_x_continuous(breaks = pretty_breaks(5)) +
scale_y_continuous(breaks = pretty_breaks(5), expand = yexpand)

if(length(data\$S_ample) != 0) {
    p <- p + facet_wrap(.~S_ample)
}

RSC
        }
        else {
            $rscript .= <<RSC;
p <- p + geom_density(alpha = 0.5, color = "#FFFFFF", size = 0.2) +
scale_x_continuous(breaks = pretty_breaks(5)) +
scale_y_continuous(breaks = pretty_breaks(5), expand = yexpand)

RSC
        }
    }

    $rscript .= <<RSC;
if(islegend == F) {
    p <- p + theme(legend.position = 'none')
}
p <- p + labs(x = xlab, y = ylab, title = title) + mytheme

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "density";
    $self->{obj} = 1;

    return $self;

}

sub count_table {
    # in order to be compatible with 1 column list

    my $file = shift;
    my $head = shift || "T";

    my %sample = ();
    my ($sample_cnt, $data_type) = (0, 0);
    my ($nrow, $ncol, $last_ncol) = (0, 0, 0);
    open my $fh, "< $file" or ERROR("count_table", "$!");
    while (<$fh>) {
        chomp;
        my @arr = split /\t/, $_;
        if($. == 1 and $head =~ /^T$|^True$/i) {
            next;
        }
        $nrow ++;
        $ncol = scalar(@arr);

        if($last_ncol == 0) {
            $last_ncol = $ncol;
        }
        elsif($last_ncol != $ncol) {
            ERROR("count_table", "the number of column isnot equal [$.]. ");
        }

        if($ncol == 1) {
            $sample_cnt = 1;
            $data_type = 1;
        }
        else {
            $sample{$arr[0]}++;
            $data_type = 2;
        }
    }
    close $fh;
    
    if($data_type == 2) {
        $sample_cnt = scalar(keys %sample);
    }

    return ($sample_cnt, $data_type);

}

return 1;

