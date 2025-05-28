package GeneralPlot::Image::Box;
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
        'type'   => 'box',
        'method' => 'none',
        'obj'    => 0,
        'save'   => 0,
        'run'    => 0,
    };

    bless $self, $class;

    return $self;

}

#-------------------------------------------------------------------------------
#  general a dot 
#-------------------------------------------------------------------------------
sub plot {

    my ($self, %opts) = @_;

    $opts{'-method'} ||= "none";

    if(defined $opts{'-method'} and $opts{'-method'} eq 'box') {
        $self->{method} = "box";
        plot_box($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'violin') {
        $self->{method} = "violin";
        plot_violin($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'splitviolin') {
        $self->{method} = "splitviolin";
        plot_splitviolin($self, %opts);
    }
    else {
        ERROR("$self->{type}", "check the -method [$opts{'-method'}]. ");
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general the box 
#-------------------------------------------------------------------------------
sub plot_box {

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

    $self->{pic}->{format}  = $opts{'-format'}  ? $opts{'-format'}  : "matrix";  # matrix | table
    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{rowcol}  = $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";     # col | row
    $self->{pic}->{group}   = $opts{'-group'}   ? $opts{'-group'}   : "none";

    $self->{pic}->{namefix} = defined $opts{'-namefix'} ? $opts{'-namefix'} : "none"; # _fpkm ...
    $self->{pic}->{filter0} = defined $opts{'-filter0'} ? $opts{'-filter0'} : "F";    # T | F
    $self->{pic}->{scale}   = defined $opts{'-scale'}   ? $opts{'-scale'}   : "none"; # log2 | log10 
    $self->{pic}->{ymin}    = defined $opts{'-ymin'}    ? $opts{'-ymin'}    : -999999;
    $self->{pic}->{ymax}    = defined $opts{'-ymax'}    ? $opts{'-ymax'}    : 999999;

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";    # 'red,blue'

    $self->{pic}->{report}  = $opts{'-report'}  ? $opts{'-report'}  : "T";

    #$self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "Box Plot";
    $self->{pic}->{title}   = defined $opts{'-title'} ? $opts{'-title'} : "Box Plot";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{adddot}  = $opts{'-adddot'} ? $opts{'-adddot'} : "F";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 6;

    # color settings
    my ($sample_cnt, $type_cnt) = (0, 0);
    if($self->{pic}->{format} eq 'table') {
        #--------------------------------------
        # two or three columns table
        # we assume the table3 is like: 
        # id    attr1   type
        # sample1   v1  t1
        # sample1   v2  t2
        # sample2   v3  t3
        # sample2   v4  t4
        #-------------------------------------
        my @return = data_check($self->{pic}->{file}, "table3");
        $sample_cnt = $return[0];
        $type_cnt   = $return[1] if(defined $return[1]);
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        if($self->{pic}->{namefix} eq 'none') {
            my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix");
            $sample_cnt = $ncol - 1;
            if($self->{pic}->{rowcol} eq 'row') {
                $sample_cnt = $nrow - 1;
            }
        }
        else {
            my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix", $self->{pic}->{namefix});
            $sample_cnt = $ncol;
        }
        
        if($self->{pic}->{group} ne "none" and -e $self->{pic}->{group}) {
            my @return = data_check($self->{pic}->{group}, "group");
            $type_cnt = $return[1];
        }
    }
    else {
        ERROR("plot_box", "check the -format [$self->{pic}->{format}] . ");
    }

    my $color = "";
    if($sample_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    if($type_cnt > 0) {
        my ($type, $tag, $num) = ("box", "default", $type_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $type_cnt) {
            $color = join(",", map {"'$_'"} @$colors[0 .. ($type_cnt-1)]);
        }
        else {
            $color = join(",", map {"'$_'"} @$colors);
            $color = "colorRampPalette(c($color))($type_cnt)";
        }
    }
    else {
        my ($type, $tag, $num) = ("box", "default", $sample_cnt);
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

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_box_theme_default();

    my $rscript = <<RSC;
$header
$theme
source("$SOFT{'ggplot2'}")

readTable <- function (file, rowcol = "col") {
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
scale   <- "$self->{pic}->{scale}"
filter0 <- as.logical($self->{pic}->{filter0})

title   <- "$self->{pic}->{title}"
xlab    <- "$self->{pic}->{xlab}"
ylab    <- "$self->{pic}->{ylab}"
legend  <- as.logical("$self->{pic}->{legend}")
adddot  <- as.logical("$self->{pic}->{adddot}")
report  <- as.logical("$self->{pic}->{report}")

ymin    <- $self->{pic}->{ymin}
ymax    <- $self->{pic}->{ymax}

width   <- $self->{pic}->{width}
height  <- $self->{pic}->{height}
colors <- c($color)

RSC

    if($self->{pic}->{format} eq 'table') {
        $rscript .= <<RSC;
data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
if(ncol(data) == 2) {
    colnames(data) <- c("S_ample", "V_alue")
} else if(ncol(data) == 3) {
    colnames(data) <- c("S_ample", "V_alue", "C_lass")
} else {
    cat("Error: the data format is irregular, check it. \\n")
    q()
}

RSC
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        $rscript .= <<RSC;
data.raw <- readTable(file, "$self->{pic}->{rowcol}")

RSC
        if($self->{pic}->{namefix} ne 'none') {
            $rscript .= <<RSC;
getcol <- grep("$self->{pic}->{namefix}", colnames(data.raw))
if(length(getcol) > 0) {
    data.raw <- data.raw[, getcol]
    colnames(data.raw) <- gsub("$self->{pic}->{namefix}", "", colnames(data.raw))
} else {
    cat("Error: no data found, check it. \\n")
    q()
}
data <- melt(data.raw)
colnames(data) <- c("S_ample", "V_alue")

RSC
        }
        else {
            $rscript .= <<RSC;
data <- melt(data.raw)
colnames(data) <- c("S_ample", "V_alue")

RSC
        }

        if($self->{pic}->{group} ne "none" and -e $self->{pic}->{group}) {
            $rscript .= <<RSC;
group <- "$self->{pic}->{group}"
class <- read.table(group, header = F, sep = '\\t', check.names = F, quote = "")
colnames(class) <- c("names", "group")
rownames(class) <- class[, 1]
class\$group <- factor(class\$group, levels = unique(class\$group))
class <- class[as.vector(data\$S_ample), ]
data\$C_lass <- class\$group

RSC
        }
    }

    $rscript .= <<RSC;
if(as.logical(filter0)) {
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

isgroup = F
if(ncol(data) == 3) {
    groups = as.character(unique(data\$T_ype))
    isgroup = T
}

min <- min(data\$V_alue)
max <- max(data\$V_alue)
sample.list = unique(data\$S_ample)

i = 0
mat.raw <- matrix(ncol = length(sample.list), nrow = 8)
if(isgroup) {
    group.list = unique(data\$C_lass)
    mat.raw = matrix(ncol = length(sample.list) * length(group.list), nrow = 9)
    for (group in group.list) {
        for (sample in sample.list) {
            tmp = data[data\$S_ample == sample & data\$C_lass == group, 2]
            if(length(tmp) > 0) {
                if(ymin != "-999999" || ymax != "999999") {
                    tmp = tmp[(tmp >= ymin & tmp <= ymax)]
                }
                quar = as.vector(quantile(tmp, probs = c(0.25,0.5,0.75)))
                i = i + 1
                mat.raw[, i] = c(sample, group, mean(tmp), min(tmp), quar, max(tmp), quar[3]-quar[1])
            } else {
                i = i + 1
                mat.raw[, i] = c(sample, group, NA, NA, NA, NA, NA, NA, NA)
            }
        }
    }
} else {
    for (sample in sample.list) {
        tmp = data[data\$S_ample == sample, 2]
        i = i + 1
        if(length(tmp) > 0) {
            if(ymin != "-999999" || ymax != "999999") {
                tmp = tmp[(tmp >= ymin & tmp <= ymax)]
            }
            quar = as.vector(quantile(tmp, probs = c(0.25,0.5,0.75)))
            mat.raw[,i] = c(sample, mean(tmp), min(tmp), quar, max(tmp), quar[3]-quar[1])
        } else {
            mat.raw[,i] = c(sample, NA, NA, NA, NA, NA, NA, NA)
        }
    }
}
mat <- t(mat.raw)
df <- na.omit(data.frame(mat))
if(isgroup) {
    colnames(df) = c("sample", "group", "mean", "min", "lower_quartiles(Q1)", "median(Q2)", 
    "upper_quartiles(Q3)", "max", "Inter-Quartile_Range(IQR)")
} else {
    colnames(df) = c("sample", "mean", "min", "lower_quartiles(Q1)", "median(Q2)", 
    "upper_quartiles(Q3)", "max", "Inter-Quartile_Range(IQR)")
}
if(report == T) {
    write.table(df, file = paste0(outdir, "/", outname, ".quantile.xls"), sep="\\t", 
        quote = FALSE, row.names = FALSE)
}

if(isgroup) {
    colnames(df) = c("sample", "group", "mean", "min", "Q1", "Q2", "Q3", "max", "IQR")
    for(n in 3:9) {
        df[, n] = as.numeric(as.character(df[, n]))
    }
} else {
    colnames(df) = c("sample", "mean", "min", "Q1", "Q2", "Q3", "max", "IQR")
    for(n in 2:8) {
        df[, n] = as.numeric(as.character(df[, n]))
    }
}

data\$S_ample = factor(data\$S_ample, levels = unique(data\$S_ample))

fig.width <- FigWidth(sample.list)
axis.size <- AxisTextSize(sample.list, width = fig.width)
if(isgroup) {
    if(length(colors) >= length(unique(data\$C_lass))) {
        colors <- colors[1:length(unique(data\$C_lass))]
    } else {
        colors <- colorRampPalette(c(colors))(length(unique(data\$C_lass)))
    }
} else if(length(colors) != length(sample.list)) {
    if(length(colors) >= length(sample.list)) {
        colors <- colors[1:length(sample.list)]
    } else {
        colors <- colorRampPalette(c(colors))(length(sample.list))
    }
}

box.width = 1
bar.width = 0.3
col.num = length(unique(data\$S_ample))
if(col.num > 100) {
    box.width = box.width / 16
} else if(col.num > 50) {
    box.width = box.width / 8
} else if(col.num > 20) {
    box.width = box.width / 4
} else if(col.num > 10) {
    box.width = box.width / 2
} else {
    box.width = box.width / 2
}

if(isgroup) {
    attr.num <- length(unique(data[data[, 1] == sample.list[1], ]\$C_lass))
    if(attr.num > 1) {
        box.width = box.width / (attr.num * 0.8)
    }
}
bar.width = box.width * 0.3

dodge = position_dodge(width = 0.7)

p <- ggplot(data = data, aes(x = S_ample, y = V_alue))
if (isgroup) {
    if(adddot) {
        p <- p + stat_boxplot(geom = "errorbar", aes(fill = C_lass), width = bar.width,
            position = dodge) + 
        geom_boxplot(aes(fill = C_lass), width = box.width, position = dodge)
    } else {
        p <- p + stat_boxplot(geom = "errorbar", aes(fill = C_lass), width = bar.width, 
            position = dodge) + 
        geom_boxplot(aes(fill = C_lass), width = box.width, position = dodge, outlier.shape = NA)
    }
} else {
    if(adddot) {
        p <- p + stat_boxplot(geom = "errorbar", width = bar.width) + 
        geom_boxplot(aes(fill = S_ample), width = box.width, position = dodge)
    } else {
        p <- p + stat_boxplot(geom = "errorbar", width = bar.width) + 
        geom_boxplot(aes(fill = S_ample), width = box.width, outlier.shape = NA, position = dodge)
    }
}
p <- p + scale_fill_manual(values = colors)
p <- p + labs(x = xlab, y = ylab, title = title)

if(ymin != -999999 && ymax != 999999) {
    p <- p + ylim(ymin, ymax)
} else if(ymin != -999999) {
    p <- p + ylim(ymin, max)
} else if(ymax != 999999) {
    p <- p + ylim(min, ymax)
}

p <- p + mytheme + xlab_angle(sample.list, width = fig.width)
ggsave(paste0(outdir, "/", outname, ".pdf"), p, width = fig.width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "box";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  general the violin 
#-------------------------------------------------------------------------------
sub plot_violin {

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

    $self->{pic}->{format}  = $opts{'-format'}  ? $opts{'-format'}  : "matrix";  # matrix | table
    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{rowcol}  = $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";     # col | row
    $self->{pic}->{group}   = $opts{'-group'}   ? $opts{'-group'}   : "none";
    
    $self->{pic}->{namefix} = defined $opts{'-namefix'} ? $opts{'-namefix'} : "none"; # _fpkm ...
    $self->{pic}->{filter0} = defined $opts{'-filter0'} ? $opts{'-filter0'} : "F";    # T | F
    $self->{pic}->{scale}   = defined $opts{'-scale'}   ? $opts{'-scale'}   : "none"; # log2 | log10 
    $self->{pic}->{ymin}    = defined $opts{'-ymin'}    ? $opts{'-ymin'}    : -999999;
    $self->{pic}->{ymax}    = defined $opts{'-ymax'}    ? $opts{'-ymax'}    : 999999;

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";    # 'red,blue'

    $self->{pic}->{report}  = $opts{'-report'}  ? $opts{'-report'}  : "T";

    #$self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "Violin Plot";
    $self->{pic}->{title}   = defined $opts{'-title'}  ? $opts{'-title'}  : "Violin Plot";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{adddot}  = $opts{'-adddot'} ? $opts{'-adddot'} : "F";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 6;

    # color settings
    my ($sample_cnt, $type_cnt) = (0, 0);
    if($self->{pic}->{format} eq 'table') {
        #--------------------------------------
        # two or three columns table
        # we assume the table3 is like: 
        # id    attr1   type
        # sample1   v1  t1
        # sample1   v2  t2
        # sample2   v3  t3
        # sample2   v4  t4
        #-------------------------------------
        my @return = data_check($self->{pic}->{file}, "table3");
        $sample_cnt = $return[0];
        $type_cnt   = $return[1] if(defined $return[1]);
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        if($self->{pic}->{namefix} eq 'none') {
            my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix");
            $sample_cnt = $ncol - 1;
            if($self->{pic}->{rowcol} eq 'row') {
                $sample_cnt = $nrow - 1;
            }
        }
        else {
            my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix", $self->{pic}->{namefix});
            $sample_cnt = $ncol;
        }

        if($self->{pic}->{group} ne "none" and -e $self->{pic}->{group}) {
            my @return = data_check($self->{pic}->{group}, "group");
            $type_cnt = $return[1];
        }
    }
    else {
        ERROR("plot_violin", "check the -format [$self->{pic}->{format}] . ");
    }

    my $color = "";
    if($sample_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    if($type_cnt > 0) {
        my ($type, $tag, $num) = ("box", "default", $type_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $type_cnt) {
            $color = join(",", map {"'$_'"} @$colors[0 .. ($type_cnt-1)]);
        }
        else {
            $color = join(",", map {"'$_'"} @$colors);
            $color = "colorRampPalette(c($color))($type_cnt)";
        }
    }
    else {
        my ($type, $tag, $num) = ("box", "default", $sample_cnt);
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
            $color = join(",", map {"'$_'"} @$colors);
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
    elsif($self->{pic}->{collist} ne "none" and $self->{pic}->{collist} =~ /,/) {
        if($self->{pic}->{collist} =~ /,/) {  
            my @colors = split /,/, $self->{pic}->{collist};
            $color = join(",", map {"'$_'"} @colors);
        }
        else {
            $color = "'$$self->{pic}->{collist}'";
        }
    }

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_box_theme_default();

    my $rscript = <<RSC;
$header
$theme
source("$SOFT{'ggplot2'}")

readTable <- function (file, rowcol = "col") {
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
scale   <- "$self->{pic}->{scale}"
filter0 <- $self->{pic}->{filter0}

title   <- "$self->{pic}->{title}"
xlab    <- "$self->{pic}->{xlab}"
ylab    <- "$self->{pic}->{ylab}"
legend  <- as.logical("$self->{pic}->{legend}")
adddot  <- as.logical("$self->{pic}->{adddot}")
report  <- as.logical("$self->{pic}->{report}")

ymin = $self->{pic}->{ymin}
ymax = $self->{pic}->{ymax}

width   <- $self->{pic}->{width}
height  <- $self->{pic}->{height}
colors <- c($color)

RSC

    if($self->{pic}->{format} eq 'table') {
        $rscript .= <<RSC;
data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
if(ncol(data) == 2) {
    colnames(data) <- c("S_ample", "V_alue")
} else if(ncol(data) == 3) {
    colnames(data) <- c("S_ample", "V_alue", "C_lass")
} else {
    cat("Error: the data format is irregular, check it. \\n")
    q()
}

RSC
    }
    elsif($self->{pic}->{format} eq 'matrix') {
        $rscript .= <<RSC;
data.raw <- readTable(file, "$self->{pic}->{rowcol}")

RSC
        if($self->{pic}->{namefix} ne 'none') {
            $rscript .= <<RSC;
getcol <- grep("$self->{pic}->{namefix}", colnames(data.raw))
if(length(getcol) > 0) {
    data.raw <- data.raw[, getcol]
    colnames(data.raw) <- gsub("$self->{pic}->{namefix}", "", colnames(data.raw))
} else {
    cat("Error: no data found, check it. \\n")
    q()
}
data <- melt(data.raw)
colnames(data) <- c("S_ample", "V_alue")

RSC
        }
        else {
            $rscript .= <<RSC;
data <- melt(data.raw)
colnames(data) <- c("S_ample", "V_alue")

RSC
        }

        if($self->{pic}->{group} ne "none" and -e $self->{pic}->{group}) {
            $rscript .= <<RSC;
group <- "$self->{pic}->{group}"
class <- read.table(group, header = F, sep = '\\t', check.names = F, quote = "")
colnames(class) <- c("names", "group")
rownames(class) <- class[, 1]
class\$group <- factor(class\$group, levels = unique(class\$group))
class <- class[as.vector(data\$S_ample), ]
data\$C_lass <- class\$group

RSC
        }
    }

    $rscript .= <<RSC;
if(as.logical(filter0)) {
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

isgroup = F
if(ncol(data) == 3) {
    groups = as.character(unique(data\$T_ype))
    isgroup = T
}

min <- min(data\$V_alue)
max <- max(data\$V_alue)
sample.list = unique(data\$S_ample)

i = 0
mat.raw <- matrix(ncol = length(sample.list), nrow = 8)
if(isgroup) {
    group.list = unique(data\$C_lass)
    mat.raw = matrix(ncol = length(sample.list) * length(group.list), nrow = 9)
    for (group in group.list) {
        for (sample in sample.list) {
            tmp = data[data\$S_ample == sample & data\$C_lass == group, 2]
            if(length(tmp) > 0) {
                if(ymin != "-999999" || ymax != "999999") {
                    tmp = tmp[(tmp >= ymin & tmp <= ymax)]
                }
                quar = as.vector(quantile(tmp, probs = c(0.25,0.5,0.75)))
                i = i + 1
                mat.raw[, i] = c(sample, group, mean(tmp), min(tmp), quar, max(tmp), quar[3]-quar[1])
            } else {
                i = i + 1
                mat.raw[, i] = c(sample, group, NA, NA, NA, NA, NA, NA, NA)
            }
        }
    }
} else {
    for (sample in sample.list) {
        tmp = data[data\$S_ample == sample, 2]
        i = i + 1
        if(length(tmp) > 0) {
            if(ymin != "-999999" || ymax != "999999") {
                tmp = tmp[(tmp >= ymin & tmp <= ymax)]
            }
            quar = as.vector(quantile(tmp, probs = c(0.25,0.5,0.75)))
            mat.raw[,i] = c(sample, mean(tmp), min(tmp), quar, max(tmp), quar[3]-quar[1])
        } else {
            mat.raw[,i] = c(sample, NA, NA, NA, NA, NA, NA, NA)
        }
    }
}
mat <- t(mat.raw)
df <- na.omit(data.frame(mat))
if(isgroup) {
    colnames(df) = c("sample", "group", "mean", "min", "lower_quartiles(Q1)", "median(Q2)", 
    "upper_quartiles(Q3)", "max", "Inter-Quartile_Range(IQR)")
} else {
    colnames(df) = c("sample", "mean", "min", "lower_quartiles(Q1)", "median(Q2)", 
    "upper_quartiles(Q3)", "max", "Inter-Quartile_Range(IQR)")
}
if(report == T) {
    write.table(df, file = paste0(outdir, "/", outname, ".quantile.xls"), sep="\\t", 
        quote = FALSE, row.names = FALSE)
}

if(isgroup) {
    colnames(df) = c("sample", "group", "mean", "min", "lower", "median", "upper", "max", "IQR")
    for(n in 3:9) {
        df[, n] = as.numeric(as.character(df[, n]))
    }
} else {
    colnames(df) = c("sample", "mean", "min", "lower", "median", "upper", "max", "IQR")
    for(n in 2:8) {
        df[, n] = as.numeric(as.character(df[, n]))
    }
}

data\$S_ample = factor(data\$S_ample, levels = unique(data\$S_ample))

fig.width <- FigWidth(sample.list)
axis.size <- AxisTextSize(sample.list, width = fig.width)
if(isgroup) {
    if(length(colors) >= length(unique(data\$C_lass))) {
        colors <- colors[1:length(unique(data\$C_lass))]
    } else {
        colors <- colorRampPalette(c(colors))(length(unique(data\$C_lass)))
    }
} else if(length(colors) != length(sample.list)) {
    if(length(colors) >= length(sample.list)) {
        colors <- colors[1:length(sample.list)]
    } else {
        colors <- colorRampPalette(c(colors))(length(sample.list))
    }
}

grouptype = 1
if(isgroup) {
    attr.num <- length(unique(data[data[, 1] == sample.list[1], ]\$C_lass))
    if(attr.num > 1) {
        grouptype = 2
    }
}

box.width = 0.1
psize = 4
col.num = length(unique(data\$S_ample))
if(grouptype == 2) col.num = col.num * length(unique(data[data[, 1] == sample.list[1], ]\$C_lass))

if(col.num > 100) {
    psize = 0.1
    box.width = box.width / 16
} else if(col.num > 50) {
    psize = 0.5
    box.width = box.width / 8
} else if(col.num > 20) {
    psize = 1
    box.width = box.width / 4
} else if(col.num > 10) {
    psize = 2
    box.width = box.width / 2
} else {
    box.width = box.width / 1.5
    psize = 2.5
}

dodge = position_dodge(width = 0.7)
p <- ggplot(data = data, aes(x = S_ample, y = V_alue)) 
if(isgroup) {
    if(grouptype == 1) {
        p <- p + geom_violin(aes(fill = C_lass), position = dodge, color = '#FFFFFF', alpha = 0.8) + 
        scale_fill_manual(values = colors) + 
        geom_boxplot(width = box.width, aes(fill = C_lass), fill = '#000000', color = '#000000',
            alpha = 0.8, outlier.shape = NA, position = dodge) + 
        geom_point(data = df, aes(x = sample, y = median), size = psize, color = "white")
    } else if(grouptype == 2) {
        p <- p + geom_violin(aes(fill = C_lass), position = dodge, color = '#FFFFFF', alpha = 0.8, 
            width = 0.7) + 
        scale_fill_manual(values = colors) + 
        geom_boxplot(width = box.width, aes(shape = C_lass), alpha = 0.8, 
            outlier.shape = NA, position = dodge, fill = '#000000', color = '#000000') + 
        stat_summary(data = data, aes(x = S_ample, y = V_alue, color = C_lass), 
            fun.y = median, geom = "point", shape = 19, size = psize * 0.8, position = dodge) +
        scale_color_manual(values = rep("#FFFFFF", length(unique(data\$C_lass))))
    }
} else {
    p <- p + geom_violin(aes(fill = S_ample), color = '#FFFFFF', alpha = 0.8) + 
    scale_fill_manual(values = colors) + 
    geom_boxplot(width = box.width, fill = "#000000", color = '#000000', alpha = 0.8, outlier.shape = NA) + 
    geom_point(data = df, aes(x = sample, y = median), size = psize, color = "white") + 
    theme(legend.position = 'none')
}

if(adddot) {
    p <- p + geom_jitter(color = "#00000066", size = 0.5)
}

p <- p + labs(x = xlab, y = ylab, title = title)

if(ymin != -999999 && ymax != 999999) {
    p <- p + ylim(ymin, ymax)
} else if(ymin != -999999) {
    p <- p + ylim(ymin, max)
} else if(ymax != 999999) {
    p <- p + ylim(min, ymax)
}

p <- p + mytheme + 
scale_x_discrete(expand = c(0.03, 0.03)) + 
xlab_angle(sample.list, width = fig.width)
ggsave(paste0(outdir, "/", outname, ".pdf"), p, width = fig.width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "violin";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  general the splitviolin 
#-------------------------------------------------------------------------------
sub plot_splitviolin {

    # can only be used on 2 attrs data

    my ($self, %opts) = @_;

    if(defined $opts{'-outprefix'}) {
        $opts{'-outdir'}  ||= dirname($opts{'-outprefix'});
        $opts{'-outname'} ||= basename($opts{'-outprefix'});
    }
    
    (!defined $self->{pic}->{file} and !defined $opts{'-file'} and 
    (!defined $opts{'-file1'} and !defined $opts{'-file2'})) and 
    ERROR("$self->{type}", "<-file>|<-file1,-file2> is not defined. ");
    !defined $self->{pic}->{outdir} and !defined $opts{'-outdir'} and
    ERROR("$self->{type}", "<-outdir> is not defined. ");
    !defined $self->{pic}->{outname} and !defined $opts{'-outname'} and
    ERROR("$self->{type}", "<-outname> is not defined. ");

    defined $opts{'-file'}    and $self->{pic}->{file}    = rel2abs($opts{'-file'});
    defined $opts{'-file1'}   and $self->{pic}->{file1}   = rel2abs($opts{'-file1'});
    defined $opts{'-file2'}   and $self->{pic}->{file2}   = rel2abs($opts{'-file2'});
    defined $opts{'-outdir'}  and $self->{pic}->{outdir}  = rel2abs($opts{'-outdir'});
    defined $opts{'-outname'} and $self->{pic}->{outname} = $opts{'-outname'};
    
    $self->{pic}->{file}  ||= "";
    $self->{pic}->{file1} ||= "";
    $self->{pic}->{file2} ||= "";

    $self->{pic}->{name1}   = defined $opts{'-name1'}   ? $opts{'-name1'}   : "group1";
    $self->{pic}->{name2}   = defined $opts{'-name2'}   ? $opts{'-name2'}   : "group2";

    $self->{pic}->{header}  = defined $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{rowcol}  = defined $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";   # col | row
    $self->{pic}->{namefix} = defined $opts{'-namefix'} ? $opts{'-namefix'} : "none";  # _fpkm ...
    $self->{pic}->{filter0} = defined $opts{'-filter0'} ? $opts{'-filter0'} : "F";     # T|F
    $self->{pic}->{scale}   = defined $opts{'-scale'}   ? $opts{'-scale'}   : "none";  # log2|log10 
    $self->{pic}->{ymin}    = defined $opts{'-ymin'}    ? $opts{'-ymin'}    : -999999;
    $self->{pic}->{ymax}    = defined $opts{'-ymax'}    ? $opts{'-ymax'}    : 999999;
    $self->{pic}->{colset}  = defined $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = defined $opts{'-collist'} ? $opts{'-collist'} : "none";  # 'red,blue'
    $self->{pic}->{col1}    = defined $opts{'-col1'}    ? $opts{'-col1'}    : "none";
    $self->{pic}->{col2}    = defined $opts{'-col2'}    ? $opts{'-col2'}    : "none";
    $self->{pic}->{alpha}   = defined $opts{'-alpha'}   ? $opts{'-alpha'}   : 0.75;
    $self->{pic}->{bartype} = defined $opts{'-bartype'} ? $opts{'-bartype'} : "split"; # none|center
    
    $self->{pic}->{report}  = defined $opts{'-report'}  ? $opts{'-report'}  : "T";

    $self->{pic}->{title}   = defined $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = defined $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = defined $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{adddot}  = defined $opts{'-adddot'} ? $opts{'-adddot'} : "F";
    $self->{pic}->{legend}  = defined $opts{'-legend'} ? $opts{'-legend'} : "T";

    $self->{pic}->{width}   = defined $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = defined $opts{'-height'} ? $opts{'-height'} : 7;

    # color settings
    my ($sample_cnt, $type_cnt) = (0, 0);
    if(-s $self->{pic}->{file}) {
        #--------------------------------------
        # two or three columns table
        # we assume the table3 is like: 
        # id    attr1   type
        # sample1   v1  t1
        # sample1   v2  t2
        # sample2   v3  t3
        # sample2   v4  t4
        #-------------------------------------
        $self->{pic}->{format} = 'table';
        my @return = data_check($self->{pic}->{file}, "table3");
        $sample_cnt = $return[0];
        $type_cnt   = $return[1] if(defined $return[1]);
    }
    elsif(-s $self->{pic}->{file1} and -s $self->{pic}->{file2}) {
        # only two matrix file could be used by this
        $self->{pic}->{format} = "matrix";
        data_check($self->{pic}->{file1}, "splitviolin");
        data_check($self->{pic}->{file2}, "splitviolin");
        my ($nrow1, $ncol1) = data_check($self->{pic}->{file1}, "matrix");
        my ($nrow2, $ncol2) = data_check($self->{pic}->{file2}, "matrix");
        if($ncol1 != $ncol2) {
            ERROR("plot_splitviolin", "the column num of the two input files are not equal.");
        }
        $sample_cnt = $ncol1 - 1;
        $type_cnt = 2;
    }
    else {
        ERROR("plot_splitviolin", "check the -format [$self->{pic}->{format}] . ");
    }

    my $color = "";
    if($sample_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    if($type_cnt == 2) {
        my ($type, $tag, $num) = ("box", "split", '0');
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        $color = join(",", map {"'$_'"} @$colors);
    }
    else {
        ERROR("$self->{type}", "the attrs for splitviolin could only be 2, check the data. ");
    }
    if($self->{pic}->{collist} ne "none" and $self->{pic}->{collist} =~ /,/) {
        my @colors = split /,/, $self->{pic}->{collist};
        $color = join(",", map {"'$_'"} @colors);
    }
    if($self->{pic}->{col1} ne "none" and $self->{pic}->{col2} ne "none") {
        my @colors = ($self->{pic}->{col1}, $self->{pic}->{col2});
        $color = join(",", map {"'$_'"} @colors);
    }

    # theme settings
    my $header = ggplot2_head();
    my $rscript = <<RSC;
$header
source("$SOFT{'all_theme'}", chdir = T)
mytheme <- box_theme_splitviolin()

source("$SOFT{"split_violin"}")
source("$SOFT{'ggplot2'}")

readTable <- function (file, header, rowcol = "col") {
    data <- read.table(file, head = header, sep = "\\t", check.names = F, 
    row.names = 1, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

file    <- "$self->{pic}->{file}"
file1   <- "$self->{pic}->{file1}"
file2   <- "$self->{pic}->{file2}"
name1   <- "$self->{pic}->{name1}"
name2   <- "$self->{pic}->{name2}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

filter0 <- as.logical("$self->{pic}->{filter0}")
header  <- as.logical("$self->{pic}->{header}")
rowcol  <- "$self->{pic}->{rowcol}"
namefix <- "$self->{pic}->{namefix}"
scale   <- "$self->{pic}->{scale}"

title   <- "$self->{pic}->{title}"
xlab    <- "$self->{pic}->{xlab}"
ylab    <- "$self->{pic}->{ylab}"
bartype <- "$self->{pic}->{bartype}"
adddot  <- as.logical("$self->{pic}->{adddot}")
legend  <- as.logical("$self->{pic}->{legend}")
report  <- as.logical("$self->{pic}->{report}")

width   <- $self->{pic}->{width}
height  <- $self->{pic}->{height}
colors  <- c($color)
alpha   <- $self->{pic}->{alpha}

ymin = $self->{pic}->{ymin}
ymax = $self->{pic}->{ymax}

data <- NULL
if(file != "") {
    data <- read.table(file, head = header, sep = "\\t", check.names = F)
    if(ncol(data) == 3) {
        colnames(data) <- c("S_ample", "V_alue", "C_lass")
    } else {
        stop("[Error]: the data format is irregular, check it. \\n")
    }

} else if(file1 != "" && file2 != "") {
    data1 <- readTable(file1, header, rowcol)
    if(namefix != "none") {
        getcol <- grep(namefix, colnames(data1))
        if(length(getcol) > 0) {
            data1 <- data1[, getcol]
            colnames(data1) <- gsub(namefix, "", colnames(data1))
        } else {
            stop("[Error]: no data found in data1, check it. \\n")
        }
    }
    data1.m <- data1 %>% melt() %>% mutate(class = name1)
    colnames(data1.m) <- c("S_ample", "V_alue", "C_lass")

    data2 <- readTable(file2, header, rowcol)
    if(namefix != "none") {
        getcol <- grep(namefix, colnames(data2))
        if(length(getcol) > 0) {
            data2 <- data2[, getcol]
            colnames(data2) <- gsub(namefix, "", colnames(data2))
        } else {
            stop("[Error]: no data found in data2, check it. \\n")
        }
    }
    data2.m <- data2 %>% melt() %>% mutate(class = name2)
    colnames(data2.m) <- c("S_ample", "V_alue", "C_lass")
    data <- rbind(data1.m, data2.m)

} else {
    stop("[Error]: check input file \\n")
}

if(filter0) {
    data <- data[data\$V_alue > 0, ]
}

# format the data
if(scale == "log2+1") {
    data\$V_alue <- log2(data\$V_alue + 1)
} else if(scale == "log10+1") {
    data\$V_alue <- log10(data\$V_alue + 1)
} else if(scale == "log2") {
    data\$V_alue[data\$V_alue == 0] = 0.0001
    data\$V_alue <- log2(data\$V_alue)
} else if(scale == "log10") {
    data\$V_alue[data\$V_alue == 0] = 0.0001
    data\$V_alue <- log10(data\$V_alue)
}

min <- min(data\$V_alue)
max <- max(data\$V_alue)
sample.list = unique(data\$S_ample)
group.list = unique(data\$C_lass)

i = 0
mat.raw = matrix(ncol = length(sample.list) * length(group.list), nrow = 9)
for (group in group.list) {
    for (sample in sample.list) {
        tmp = data[data\$S_ample == sample & data\$C_lass == group, 2]
        if(length(tmp) > 0) {
            if(ymin != "-999999" || ymax != "999999") {
                tmp = tmp[(tmp >= ymin & tmp <= ymax)]
            }
            quar = as.vector(quantile(tmp, probs = c(0.25,0.5,0.75)))
            i = i + 1
            mat.raw[, i] = c(sample, group, mean(tmp), min(tmp), quar, max(tmp), quar[3]-quar[1])
        } else {
            i = i + 1
            mat.raw[, i] = c(sample, group, NA, NA, NA, NA, NA, NA, NA)
        }
    }
}

mat <- t(mat.raw)
df <- na.omit(data.frame(mat))

colnames(df) = c("sample", "group", "mean", "min", "lower_quartiles(Q1)", "median(Q2)", 
"upper_quartiles(Q3)", "max", "Inter-Quartile_Range(IQR)")
if(report == T) {
    write.table(df, file = paste0(outdir, "/", outname, ".quantile.xls"), sep="\\t", 
        quote = FALSE, row.names = FALSE)
}

colnames(df) = c("sample", "group", "mean", "min", "lower", "median", "upper", "max", "IQR")
for(n in 3:9) {
    df[, n] = as.numeric(as.character(df[, n]))
}
data\$S_ample = factor(data\$S_ample, levels = unique(data\$S_ample))

fig.width <- FigWidth(sample.list)
axis.size <- AxisTextSize(sample.list, width = fig.width)
box.width <- 0.1
psize <- 2
col.num <- length(unique(data\$S_ample))
if(col.num > 100) {
    box.width = box.width / 16
    psize = psize / 10
} else if(col.num > 50) {
    box.width = box.width / 8
    psize = psize / 6
} else if(col.num > 20) {
    box.width = box.width / 4
    psize = psize / 4
} else if(col.num > 10) {
    box.width = box.width / 2
    psize = psize / 2
} else {
    box.width = box.width
    psize = psize
}

p <- ggplot(data = data, aes(x = S_ample, y = V_alue, fill = C_lass)) + 
geom_split_violin(trim = FALSE, color = "#FFFFFF", alpha = alpha)
if(bartype == "center" || bartype == "centre") {
    p <- p + geom_boxplot(width = box.width, aes(ymin = ..lower.., ymax = ..upper.., shape = C_lass),
            outlier.shape = NA, position = position_dodge(box.width),
            fill = "#000000", color = '#FFFFFF', size = 0.1, alpha = 0.8) +
        geom_point(data = df, aes(x = sample, y = median, fill = group),
            position = position_dodge2(box.width), color = "#FFFFFF", size = psize)
} else if(bartype == "split") {
    p <- p + geom_boxplot(width = box.width, aes(ymin = ..lower.., ymax = ..upper.., shape = C_lass), 
            outlier.shape = NA, position = position_dodge(0.5), 
            fill = "#000000", color = NA, size = 0.1, alpha = 0.8) + 
        geom_point(data = df, aes(x = sample, y = median, fill = group), 
            position = position_dodge2(0.5), color = "#FFFFFF", size = psize)
}
if(adddot) {
    p <- p + geom_jitter(color = "#00000099", size = 0.5)
}
p <- p + scale_fill_manual(values = colors) + 
labs(x = xlab, y = ylab, title = title)

if(ymin != -999999 && ymax != 999999) {
    p <- p + ylim(ymin, ymax)
} else if(ymin != -999999) {
    p <- p + ylim(ymin, max)
} else if(ymax != 999999) {
    p <- p + ylim(min, ymax)
}

p <- p + mytheme + 
xlab_angle(sample.list, width = fig.width)
ggsave(paste0(outdir, "/", outname, ".pdf"), p, width = fig.width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "splitviolin";
    $self->{obj} = 1;

    return $self;

}


return 1;

