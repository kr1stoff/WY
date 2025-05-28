package GeneralPlot::Image::Pie;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.10 16:11:07          |
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
        'type'   => 'pie',
        'method' => 'none',
        'obj'    => 0,
        'save'   => 0,
        'run'    => 0,
    };

    bless $self, $class;

    return $self;

}

#-------------------------------------------------------------------------------
#  general a single pie chart from a format table directly
#-------------------------------------------------------------------------------
sub plot {

    my ($self, %opts) = @_;

    if(defined $opts{'-method'} and $opts{'-method'} eq 'twopie') {
        $self->{method} = "twopie";
        plot_twopie($self, %opts);
    }
    else {
        $self->{method} = "pie";
        plot_pie($self, %opts);
    }

    return $self;
}

sub plot_pie {

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
    $self->{pic}->{colrev}  = $opts{'-colrev'}  ? $opts{'-colrev'}  : "T";

    $self->{pic}->{order}     = $opts{'-order'}     ? $opts{'-order'}     : "none"; # none | ab | ba
    $self->{pic}->{direction} = $opts{'-direction'} ? $opts{'-direction'} : -1; # 1 | -1 
    $self->{pic}->{header}  = $opts{'-header'} ? $opts{'-header'} : "T";
    $self->{pic}->{label}   = $opts{'-label'}  ? $opts{'-label'}  : "F";
    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";
    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 7;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 7;

    # color
    my $color = "";
    my @return = data_check($self->{pic}->{file}, "table3");
    my $label_cnt = $return[0];
    $label_cnt ++ if($self->{pic}->{header} !~ /^T$|^True$/i);
    if($label_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("pie", "default", $label_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }
        
        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if($self->{pic}->{colrev} =~ /^T$|^True$/) {
            @$colors = reverse(@$colors);
        }
        if(scalar(@$colors) >= $label_cnt) {
            $color = join(",", map {"'$_'"} @$colors[0 .. ($label_cnt-1)]);
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
    elsif($self->{pic}->{collist} ne 'none') {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            if(scalar(@colors) != $label_cnt) {
                my $col_cnt = scalar(@colors);
                warn("Insufficient colors. $label_cnt needed but only $col_cnt provided. \n");
            }
            $color = join(",", map {"'$_'"} @colors);
        }
    }

    # theme
    my $header = ggplot2_head();
    my $theme  = ggplot2_pie_theme_default();

    my $rscript = <<RSC;
$header
$theme

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

title  <- "$self->{pic}->{title}"
xlab   <- "$self->{pic}->{xlab}"
ylab   <- "$self->{pic}->{ylab}"
islegend <- as.logical("$self->{pic}->{legend}")
islabel  <- as.logical("$self->{pic}->{label}")

width = $self->{pic}->{width}
height = $self->{pic}->{height}
colors <- c($color)

order <- as.character("$self->{pic}->{order}")
direction <- as.character("$self->{pic}->{direction}")

if(direction == "clockwise" || direction == "1") {
    direction = 1
} else if(direction == "anticlockwise" || direction == "-1") {
    direction = -1
} else {
    direction = -1
}

data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", quote = "", check.names = F)
if(length(data) == 2) {
    colnames(data) <- c("x_index", "y_count")
    tag.name <- data\$x_index
    #data <- data %>% arrange(desc(y_count))  # cannot sort because of uniformity in samples
    if(order == "ab") {
        data <- data %>% arrange(y_count, x_index)
    } else if(order == "ba") {
        data <- data %>% arrange(desc(y_count), x_index)
    }
    data\$x_index <- factor(data\$x_index, levels = unique(data\$x_index), ordered = T)

    data.melt <- melt(data, id.vars = c("x_index"))
    data.melt.sum <- sum(data.melt\$value)
    data.melt\$P_ercent <- data.melt\$value / data.melt.sum
    data.melt\$P_ercent2 <- sprintf("%0.1f%%", data.melt\$P_ercent * 100)
    data.melt\$L_abel <- sprintf("%s (%0.2f%%)", data.melt\$x_index, data.melt\$P_ercent * 100)
    data.melt\$L_abel <- factor(data.melt\$L_abel, levels = unique(data.melt\$L_abel), ordered = T)

    p <- ggplot(data.melt, aes(x = "", y = P_ercent, fill = x_index, color = I("#FFFFFF"))) +
    geom_bar(stat = "identity", width = 1, size = 0.2, position = position_stack())

    if(islabel == T) {
        p <- p + geom_text(aes(x = 1.25, label = P_ercent2), color = "#000000", size = 2.5,
        position = position_stack(reverse = F, vjust = 0.5), hjust = 0.5)
    }

    label.num <- length(unique(data.melt\$x_index))
    if(length(colors) >= label.num) {
        colors <- colors[1:label.num]
    } else {
        colors <- colorRampPalette(colors)(label.num)
    }
    if(order == "ba") {
        colors = rev(colors)
    }
    p <- p + scale_fill_manual(values = colors, labels = data.melt\$L_abel) +
    coord_polar(theta = "y", direction = direction) +
    labs(x = "", y = "", title = title) + 
    mytheme + 
    theme(plot.title = element_text(vjust = -2))

} else if(length(data) == 3) {
    colnames(data)[1] = "x_index"
    colnames(data)[2] = "y_count"
    colnames(data)[3] = "z_object"
    if(order == "ab") {
        data <- data %>% arrange(z_object, y_count)
    } else if(order == "ba") {
        data <- data %>% arrange(z_object, desc(y_count))
    }
    data\$x_index <- factor(data\$x_index, levels = unique(data\$x_index), ordered = T)
    data\$z_object <- factor(data\$z_object, levels = unique(data\$z_object), ordered = T)

    data.melt <- melt(data, id.vars = c("x_index", "z_object"))
    data.melt.sum <- aggregate(data.melt\$value, list(z_object = data.melt\$z_object), sum)
    rownames(data.melt.sum) <- data.melt.sum\$z_object
    data.melt\$P_ercent <- data.melt\$value / data.melt.sum[data.melt\$z_object, "x"]
    data.melt\$P_ercent2 <- sprintf("%0.1f%%", data.melt\$P_ercent * 100)
    data.melt\$L_abel <- sprintf("%s (%0.2f%%)", data.melt\$x_index, data.melt\$P_ercent * 100)
    data.melt\$L_abel <- factor(data.melt\$L_abel, levels = unique(data.melt\$L_abel), ordered = T)

    p <- ggplot(data.melt, aes(x = "", y = P_ercent, fill = x_index, color = I("#FFFFFF"))) +
    geom_bar(stat = "identity", width = 1, size = 0.2, position = position_stack())

    if(islabel == T) {
        p <- p + geom_text(aes(x = 1.25, label = P_ercent2), color = "#000000", size = 2.8,
        position = position_stack(reverse = F, vjust = 0.5), hjust = 0.5)
    }

    label.num <- length(unique(data.melt\$x_index))
    if(length(colors) >= label.num) {
        colors <- colors[1:label.num]
    } else {
        colors <- colorRampPalette(colors)(label.num)
    }
    if(order == "ba") {
        colors = rev(colors)
    }
    p <- p + scale_fill_manual(values = colors) +
    coord_polar(theta = "y", direction = direction) +
    facet_wrap(~z_object) +
    labs(x = "", y = "", title = title) +
    mytheme +
    theme(legend.position = "bottom", plot.title = element_text(vjust = 1))
}

if(islegend != T) {
    p <- p + theme(legend.position = "none")
}

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "pie";
    $self->{obj} = 1;

    return $self;

}

sub plot_twopie {

    my ($self, %opts) = @_;

    if(defined $opts{'-outprefix'}) {
        $opts{'-outdir'}  ||= dirname($opts{'-outprefix'});
        $opts{'-outname'} ||= basename($opts{'-outprefix'});
    }
    !defined $self->{pic}->{file} and !defined $opts{'-file'} and
    ERROR("line", "<-file> is not defined. ");
    !defined $self->{pic}->{outdir} and !defined $opts{'-outdir'} and
    ERROR("line", "<-outdir> is not defined. ");
    !defined $self->{pic}->{outname} and !defined $opts{'-outname'} and
    ERROR("line", "<-outname> is not defined. ");

    defined $opts{'-file'}    and $self->{pic}->{file}    = rel2abs($opts{'-file'});
    defined $opts{'-outdir'}  and $self->{pic}->{outdir}  = rel2abs($opts{'-outdir'});
    defined $opts{'-outname'} and $self->{pic}->{outname} = $opts{'-outname'};

    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";

    $self->{pic}->{header}  = $opts{'-header'} ? $opts{'-header'} : "T";
    $self->{pic}->{label}   = $opts{'-label'}  ? $opts{'-label'}  : "F";
    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";
    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 7;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 7;

    # color
    my $color = "";
    my @return = data_check($self->{pic}->{file}, "table3");
    my ($label_cnt, $class_cnt) = @return;
    #my $color_cnt = $label_cnt + $class_cnt;
    my $color_cnt = $label_cnt;
    if($color_cnt < 1) {
        ERROR("two_pie", "no data found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("twopie", "default", $label_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $color_cnt) {
            $color = join(",", map {"'$_'"} @$colors[0 .. ($color_cnt-1)]);
        }
        else {
            $color = join(",", map { "'$_'" } @$colors);
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
            if(scalar(@colors) != $color_cnt) {
                my $col_cnt = scalar(@colors);
                warn("Insufficient colors. $color_cnt needed but only $col_cnt provided. \n");
            }
            $color = join(",", map {"'$_'"} @colors);
        }
    }

    # theme
    my $header = ggplot2_head();
    my $theme  = ggplot2_pie_theme_default();

    my $rscript = <<RSC;
$header
$theme

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

title  <- "$self->{pic}->{title}"
xlab   <- "$self->{pic}->{xlab}"
ylab   <- "$self->{pic}->{ylab}"
islabel  <- as.logical("$self->{pic}->{label}")
islegend <- as.logical("$self->{pic}->{legend}")

width = $self->{pic}->{width}
height = $self->{pic}->{height}
colors <- c($color)

data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F)
colnames(data) <- c("x_index", 'y_count', 'z_type')
label1.cnt <- length(unique(data\$z_type))
label2.cnt <- length(unique(data\$x_index))

color1 <- colorRampPalette(colors)(label1.cnt)
color1 <- as.vector(sapply(color1, function(str){paste0(str, "99")}))
color2 <- colorRampPalette(colors)(label2.cnt)
colors <- c(color1, color2)

data.sum <- data %>% group_by(z_type) %>% summarise(group.sum = sum(y_count)) %>% 
arrange(desc(group.sum))
data\$z_type <- factor(data\$z_type, levels = data.sum\$z_type, order = T)
data <- data %>% arrange(z_type, desc(y_count), x_index)

data.merge <- rbind(
    data.frame(type = "A", x = data.sum\$z_type, y = data.sum\$group.sum),
    data.frame(type = "B", x = data\$x_index, y = data\$y_count)
)
data.merge\$type <- factor(data.merge\$type, levels = unique(data.merge\$type), order = T)
data.merge\$x <- factor(data.merge\$x, levels = unique(data.merge\$x), order = T)
data.merge\$percent <- data.merge\$y / sum(data\$y_count)
data.merge\$label <- sprintf("%s (%0.2f%%)", data.merge\$x, data.merge\$percent * 100)

p <- ggplot(data.merge, aes(x = type, y = y, fill = x, color = I("#FFFFFF"))) +
geom_bar(stat = "identity", width = 1, size = 0.2, position = position_stack()) +
scale_fill_manual(values = colors, labels = data.merge\$label) +
coord_polar(theta = "y", direction = -1) +
labs(x = "", y = "", title = title) + 
mytheme

if(islabel == T) {
    p <- p + geom_text(aes(x = type, y = y, label = x), color = "#000000", size = 3,
        position = position_stack(reverse = F, vjust = 0.5), hjust = 0.5)
}

if(islegend != T) {
    p <- p + theme(legend.position = "none")
}

ggsave(file = paste0(outdir, "/", outname, ".pdf"), plot = p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "twopie";
    $self->{obj} = 1;

    return $self;

}

return 1;

