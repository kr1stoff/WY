package GeneralPlot::Image::Dot;
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
        'type'   => 'dot',
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

    if(defined $opts{'-method'} and $opts{'-method'} eq 'pca') {
        $self->{method} = "pca";
        plot_pca($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'volcano') {
        $self->{method} = "volcano";
        plot_volcano($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'manhattan') {
        $self->{method} = "manhattan";
        plot_manhattan($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'dot') {
        $self->{method} = "dot";
        plot_dot($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'scatter3d') {
        $self->{method} = 'scatter3d';
        plot_scatter3d($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'scatter_bubble') {
        $self->{method} = 'scatter_bubble';
        plot_scatter_bubble($self, %opts);
    }
    elsif(defined $opts{'-method'} and $opts{'-method'} eq 'lollipop') {
        $self->{method} = 'lollipop';
        plot_lollipop($self, %opts);
    }
    else {
        ERROR("heatmap", "check the -method [$opts{'-method'}]. ");
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general the dot 
#-------------------------------------------------------------------------------
sub plot_dot {

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
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";   # 'red,blue'

    $self->{pic}->{x}       = $opts{'-x'}       ? $opts{'-x'}      : 1;
    $self->{pic}->{y}       = $opts{'-y'}       ? $opts{'-y'}      : 2;
    $self->{pic}->{group}   = $opts{'-group'}   ? $opts{'-group'}  : 'none';
    
    $self->{pic}->{header}  = $opts{'-header'} ? $opts{'-header'} : "T";
    $self->{pic}->{title}   = defined $opts{'-title'}  ? $opts{'-title'}  : "Scatter Plot";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";
    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 7;

    # color settings
    my $color = "";
    my ($type, $tag, $num) = ("dot", "default", 10);
    if($self->{pic}->{colset} ne 'none') {
        my @set = split /\_/, $self->{pic}->{colset};
        $type = $set[0];
        $tag  = $set[1];
        $num  = defined $set[2] ? $set[2] : $num;
    }

    my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
    $color = join(",", map {"'$_'"} @$colors);

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
            $color = join(",", map {"'$_'"} @colors);
        }
        else {
            $color = "'$self->{pic}->{collist}'";
        }
    }

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_dot_theme_default();

    my $rscript = <<RSC;
$header
$theme

source("$SOFT{ggplot2}")

file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

x <- $self->{pic}->{x}
y <- $self->{pic}->{y}
group <- "$self->{pic}->{group}"

title  <- "$self->{pic}->{title}"
xlab   <- "$self->{pic}->{xlab}"
ylab   <- "$self->{pic}->{ylab}"
legend <- as.logical("$self->{pic}->{legend}")

width  <- $self->{pic}->{width}
height <- $self->{pic}->{height}
colors <- c($color)

data.r <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, quote = "")
colname <- colnames(data.r)
isgroup = T
if(group != "none") {
    group <- as.numeric(group)
    if(length(colname) >= group) {
        data <- data.frame("x" = data.r[, x], "y" = data.r[, y], "group" = data.r[, group])
    } else {
        cat("Error: the group column is not exist. \\n")
        q()
    }
} else {
    data <- data.frame("x" = data.r[, x], "y" = data.r[, y], "group" = "none")
    isgroup = F
}

p <- ggplot(data, aes(x = x, y = y, color = group))

if(is.numeric(data\$group) == T) {
    p <- p + geom_point(size = 1, alpha = 0.8) + 
    scale_color_gradientn(colors = colors)
} else {
    group.num <- length(unique(data\$group))
    if(group.num <= length(colors)) {
        colors = colors[1:group.num]
    } else {
        colors = colorRampPalette(colors)(group.num)
    }
    p <- p + geom_point(size = 1, alpha = 0.8) + 
    scale_color_manual(values = colors)
}

p <- p + labs(x = xlab, y = ylab, title = title, color = colname[group]) +
mytheme + theme(legend.title = element_text(), legend.title.align = 0)

if(is.numeric(data\$x) == F) {
    axis.size <- AxisTextSize(data\$x, width = width)
    p <- p + theme(axis.text = element_text(size = axis.size)) + 
    xlab_angle(unique(data\$x), width = width)
}

if(legend == F || isgroup == F) {
    p <- p + theme(legend.position = "none")
}

outfig <- paste0(outdir, "/", outname, ".pdf")
ggsave(outfig, p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "dot";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
# plot lollipop with a table and a group file
#-------------------------------------------------------------------------------
sub plot_lollipop {

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

    $self->{pic}->{group}        = $opts{'-group'}        ? rel2abs($opts{'-group'}) : "none";
    $self->{pic}->{scale}        = $opts{'-scale'}        ? $opts{'-scale'} : "none";
    $self->{pic}->{dot_size}     = $opts{'-dot_size'}     ? $opts{'-dot_size'} : 4;
    $self->{pic}->{line_size}    = $opts{'-line_size'}    ? $opts{'-line_size'} : 1;
    $self->{pic}->{dot_color}    = $opts{'-dot_color'}    ? $opts{'-dot_color'} : "none";
    $self->{pic}->{line_color}   = $opts{'-line_color'}   ? $opts{'-line_color'} : "grey50";
    $self->{pic}->{transparency} = $opts{'-transparency'} ? $opts{'-transparency'} : 1;
    $self->{pic}->{is_label}     = $opts{'-is_label'}     ? $opts{'-is_label'} : "F";
    $self->{pic}->{is_baseline}  = $opts{'-is_baseline'}  ? $opts{'-is_baseline'} : "F";
    $self->{pic}->{is_flip}      = $opts{'-is_flip'}      ? $opts{'-is_flip'} : "F";
    $self->{pic}->{is_sort}      = $opts{'-is_sort'}      ? $opts{'-is_sort'} : "F";
    $self->{pic}->{x_lab}        = $opts{'-x_lab'}        ? $opts{'-x_lab'} : "none";
    $self->{pic}->{y_lab}        = $opts{'-y_lab'}        ? $opts{'-y_lab'} : "none";
    $self->{pic}->{title}        = $opts{'-title'}        ? $opts{'-title'} : "none";
    $self->{pic}->{lgd_lab}      = $opts{'-lgd_lab'}      ? $opts{'-lgd_lab'} : "none";

    $self->{pic}->{scale} = "'$self->{pic}->{scale}'";
    $self->{pic}->{dot_color} = "'$self->{pic}->{dot_color}'";
    $self->{pic}->{line_color} = "'$self->{pic}->{line_color}'";
    $self->{pic}->{x_lab} = "'$self->{pic}->{x_lab}'";
    $self->{pic}->{y_lab} = "'$self->{pic}->{y_lab}'";
    $self->{pic}->{lgd_lab} = "'$self->{pic}->{lgd_lab}'";
    $self->{pic}->{title} = "'$self->{pic}->{title}'";
    
    my $shell = "export LANG='zh_CN.UTF-8' ; mkdir -p $self->{pic}->{outdir} ; ".
    "$SOFT{Rscript} $SOFT{lollipop} $self->{pic}->{file} ".
    "$self->{pic}->{outdir}/$self->{pic}->{outname} $self->{pic}->{group} ".
    "$self->{pic}->{scale} $self->{pic}->{dot_size} $self->{pic}->{line_size} ".
    "$self->{pic}->{dot_color} $self->{pic}->{line_color} $self->{pic}->{transparency} ".
    "$self->{pic}->{is_label} $self->{pic}->{is_baseline} $self->{pic}->{is_flip} ".
    "$self->{pic}->{is_sort} $self->{pic}->{x_lab} $self->{pic}->{y_lab} $self->{pic}->{title} ".
    "$self->{pic}->{lgd_lab}";

    $self->{rscript} = "";
    $self->{command} = $shell;
    $self->{method} = "lollipop";
    $self->{obj} = 1;

    return $self;

}



#-------------------------------------------------------------------------------
# plot bubble with two matrix
#-------------------------------------------------------------------------------
sub plot_scatter_bubble {

    my ($self, %opts) = @_;

    if(defined $opts{'-outprefix'}) {
        $opts{'-outdir'}  ||= dirname($opts{'-outprefix'});
        $opts{'-outname'} ||= basename($opts{'-outprefix'});
    }

    (!defined $opts{'-file1'} and !defined $opts{'-file2'}) and
    ERROR("$self->{type}", "<-file>|<-file1,-file2> is not defined. ");
    !defined $self->{pic}->{outdir} and !defined $opts{'-outdir'} and
    ERROR("$self->{type}", "<-outdir> is not defined. ");
    !defined $self->{pic}->{outname} and !defined $opts{'-outname'} and
    ERROR("$self->{type}", "<-outname> is not defined. ");

    defined $opts{'-file1'}   and $self->{pic}->{file1}   = rel2abs($opts{'-file1'});
    defined $opts{'-file2'}   and $self->{pic}->{file2}   = rel2abs($opts{'-file2'});
    defined $opts{'-outdir'}  and $self->{pic}->{outdir}  = rel2abs($opts{'-outdir'});
    defined $opts{'-outname'} and $self->{pic}->{outname} = $opts{'-outname'};

    $self->{pic}->{rowname_file} = $opts{'-rowname_file'} ? rel2abs($opts{'-rowname_file'}) : "none";
    $self->{pic}->{colname_file} = $opts{'-colname_file'} ? rel2abs($opts{'-colname_file'}) : "none";
    $self->{pic}->{scale1} = $opts{'-scale1'} ? $opts{'-scale1'} : "none";  # none | log2 | log10
    $self->{pic}->{scale2} = $opts{'-scale2'} ? $opts{'-scale2'} : "none";  # none | -log2 | -log10
    $self->{pic}->{dot_minsize} = $opts{'-dot_minsize'} ? $opts{'-dot_minsize'} : 1;
    $self->{pic}->{dot_maxsize} = $opts{'-dot_maxsize'} ? $opts{'-dot_maxsize'} : 6;
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";   # 'red,blue'
    $self->{pic}->{transparency} = $opts{'-transparency'} ? $opts{'-transparency'} : 0.9;
    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "none";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "none";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "none";
    $self->{pic}->{lgd1}    = $opts{'-lgd1'}   ? $opts{'-lgd1'}   : "none";
    $self->{pic}->{lgd2}    = $opts{'-lgd2'}   ? $opts{'-lgd2'}   : "none";

    
    $self->{pic}->{scale1}  = "'$self->{pic}->{scale1}'";
    $self->{pic}->{scale2}  = "'$self->{pic}->{scale2}'";
    $self->{pic}->{collist} = "'$self->{pic}->{collist}'";
    $self->{pic}->{title}   = "'$self->{pic}->{title}'";
    $self->{pic}->{xlab}    = "'$self->{pic}->{xlab}'";
    $self->{pic}->{ylab}    = "'$self->{pic}->{ylab}'";
    $self->{pic}->{lgd1}    = "'$self->{pic}->{lgd1}'";
    $self->{pic}->{lgd2}    = "'$self->{pic}->{lgd2}'";

    my $shell = "export LANG='zh_CN.UTF-8' ; mkdir -p $self->{pic}->{outdir} ; ".
    "$SOFT{Rscript} $SOFT{scatter_bubble} $self->{pic}->{file1} $self->{pic}->{file2} ".
    "$self->{pic}->{outdir}/$self->{pic}->{outname} ".
    "$self->{pic}->{rowname_file} $self->{pic}->{colname_file} ".
    "$self->{pic}->{scale1} $self->{pic}->{scale2} $self->{pic}->{dot_minsize} ".
    "$self->{pic}->{dot_maxsize} $self->{pic}->{collist} $self->{pic}->{transparency} ".
    "$self->{pic}->{lgd1} $self->{pic}->{lgd2} $self->{pic}->{xlab} $self->{pic}->{ylab} ".
    "$self->{pic}->{title} ";

    $self->{rscript} = "";
    $self->{command} = $shell;
    $self->{method} = "scatter_bubble";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
# plot 3D scatter with a table
#-------------------------------------------------------------------------------
sub plot_scatter3d {

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
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";   # 'red,blue'
 
    $self->{pic}->{x}         = $opts{'-x'}         ? $opts{'-x'}         : 2;
    $self->{pic}->{y}         = $opts{'-y'}         ? $opts{'-y'}         : 3;
    $self->{pic}->{z}         = $opts{'-z'}         ? $opts{'-z'}         : 4;
    $self->{pic}->{samplecol} = $opts{'-samplecol'} ? $opts{'-samplecol'} : 1;
    $self->{pic}->{groupcol}  = $opts{'-groupcol'}  ? $opts{'-groupcol'}  : 'none';
    
    $self->{pic}->{header}  = $opts{'-header'} ? $opts{'-header'} : "T";
    $self->{pic}->{title}   = $opts{'-title'}  ? $opts{'-title'}  : "";
    $self->{pic}->{xlab}    = $opts{'-xlab'}   ? $opts{'-xlab'}   : "";
    $self->{pic}->{ylab}    = $opts{'-ylab'}   ? $opts{'-ylab'}   : "";
    $self->{pic}->{zlab}    = $opts{'-zlab'}   ? $opts{'-zlab'}   : "";
    $self->{pic}->{pcex}    = $opts{'-pcex'}   ? $opts{'-pcex'}   : 1;
    $self->{pic}->{tcex}    = $opts{'-tcex'}   ? $opts{'-tcex'}   : 0.8;
    $self->{pic}->{vjust}   = $opts{'-vjust'}  ? $opts{'-vjust'}  : 0.5;

    $self->{pic}->{grid}    = $opts{'-grid'}   ? $opts{'-grid'}   : 'F';
    $self->{pic}->{label}   = $opts{'-label'}  ? $opts{'-label'}  : 'F';

    # color settings
    my $color = "";
    my ($type, $tag, $num) = ("dot", "default", 10);

    my $cnt1 = 0; my %cnt2 = ();
    open my $fh, "< $self->{pic}->{file}" or die "$!";
    while (<$fh>) {
        chomp;
        $. == 1 and next;
        my @arr = split /\t/, $_;
        $cnt1 ++;
        if($self->{pic}->{groupcol} !~ /none/ and $self->{pic}->{groupcol} =~ /^\d+$/) {
            $cnt2{$arr[$self->{pic}->{groupcol} - 1]} ++;
        }
    }
    close $fh;
    $num = scalar(keys %cnt2) > 0 ? scalar(keys %cnt2) : $cnt1;

    my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
    $color = join(",", @$colors);

    if($self->{pic}->{colnpg} =~ /^T$|^True$/i) {
        my $colors = fetch_color(-type => "ggsci", -tag => 'npg', -num => 0);
        $color = join(",", @$colors);
    }
    elsif($self->{pic}->{collanc} =~ /^T$|^True$/i) {
        my $colors = fetch_color(-type => "ggsci", -tag => 'lancet', -num => 0);
        $color = join(",", @$colors);
    }
    elsif($self->{pic}->{collist} ne 'none') {
        if($self->{pic}->{collist} =~ /,/) {
            my @colors = split /,/, $self->{pic}->{collist};
            $color = join(",", @colors);
        }
        else {
            $color = "$self->{pic}->{collist}";
        }
    }

    my $title = $self->{pic}->{title} ne "" ? "--title $self->{pic}->{title}" : "";
    my $xlab  = $self->{pic}->{xlab}  ne "" ? "--xlab  $self->{pic}->{xlab}"  : "";
    my $ylab  = $self->{pic}->{ylab}  ne "" ? "--ylab  $self->{pic}->{ylab}"  : "";
    my $zlab  = $self->{pic}->{zlab}  ne "" ? "--zlab  $self->{pic}->{zlab}"  : "";
    my $group = $self->{pic}->{groupcol} !~ /none/ ? "-u $self->{pic}->{groupcol}" : "";

    my $shell = "$SOFT{Rscript} $SOFT{scatter_3d} --file $self->{pic}->{file} ".
    "--outdir $self->{pic}->{outdir} --outname $self->{pic}->{outname} ".
    "-x $self->{pic}->{x} -y $self->{pic}->{y} -z $self->{pic}->{z} ".
    "-s $self->{pic}->{samplecol} $group $title $xlab $ylab $zlab ".
    "--pcex $self->{pic}->{pcex} --tcex $self->{pic}->{tcex} --vjust $self->{pic}->{vjust} ".
    "--color '$color' --isgrid $self->{pic}->{grid} --islabel $self->{pic}->{label}";

    $self->{rscript} = "";
    $self->{command} = $shell;
    $self->{method} = "scatter3d";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  general the pca  
#-------------------------------------------------------------------------------
sub plot_pca {

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

    defined $opts{'-file'}      and $self->{pic}->{file}      = rel2abs($opts{'-file'});
    defined $opts{'-outdir'}    and $self->{pic}->{outdir}    = rel2abs($opts{'-outdir'});
    defined $opts{'-outname'}   and $self->{pic}->{outname}   = $opts{'-outname'};
    defined $opts{'-group'}     and $self->{pic}->{group}     = rel2abs($opts{'-group'});
    defined $opts{'-arrowfile'} and $self->{pic}->{arrowfile} = rel2abs($opts{'-arrowfile'});
    
    $self->{pic}->{group} ||= 'none';
    $self->{pic}->{arrowfile} ||= 'none';

    $self->{pic}->{colnpg}   = $opts{'-colnpg'}   ? $opts{'-colnpg'}   : "F";
    $self->{pic}->{collanc}  = $opts{'-collanc'}  ? $opts{'-collanc'}  : "F";
    $self->{pic}->{colset}   = $opts{'-colset'}   ? $opts{'-colset'}   : "none";
    $self->{pic}->{collist}  = $opts{'-collist'}  ? $opts{'-collist'}  : "none";   # 'red,blue'

    $self->{pic}->{skipstat} = $opts{'-skipstat'} ? $opts{'-skipstat'} : "F";
    $self->{pic}->{pc1}      = $opts{'-pc1'}      ? $opts{'-pc1'}      : 1;
    $self->{pic}->{pc2}      = $opts{'-pc2'}      ? $opts{'-pc2'}      : 2;
    
    $self->{pic}->{arrow}    = $opts{'-arrow'}    ? $opts{'-arrow'}    : "F";
    $self->{pic}->{ellipse}  = $opts{'-ellipse'}  ? $opts{'-ellipse'}  : "F";
    $self->{pic}->{ellipse}  = $opts{'-ellipse'}  ? $opts{'-ellipse'}  : "F";
    $self->{pic}->{dotsize}  = $opts{'-dotsize'}  ? $opts{'-dotsize'}  : "none";

    $self->{pic}->{report}  = $opts{'-report'}  ? $opts{'-report'}  : "T";

    $self->{pic}->{rowcol}  = $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";    # row
    $self->{pic}->{namefix} = $opts{'-namefix'} ? $opts{'-namefix'} : "none";   # _fpkm ...

    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'}  : "T";
    $self->{pic}->{label}   = $opts{'-label'}   ? $opts{'-label'}   : "F";
    $self->{pic}->{title}   = defined $opts{'-title'}   ? $opts{'-title'}   : "PCA";
    $self->{pic}->{xlab}    = $opts{'-xlab'}    ? $opts{'-xlab'}    : "none";
    $self->{pic}->{ylab}    = $opts{'-ylab'}    ? $opts{'-ylab'}    : "none";
    $self->{pic}->{width}   = $opts{'-width'}   ? $opts{'-width'}   : 8;
    $self->{pic}->{height}  = $opts{'-height'}  ? $opts{'-height'}  : 8;
    $self->{pic}->{legend}  = $opts{'-legend'}  ? $opts{'-legend'}  : "T";

    # color settings
    my $color = "";
    my $sample_cnt = 0;
    if($self->{pic}->{namefix} eq 'none') {
        my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix");
        $sample_cnt = $ncol - 1;
        if($self->{pic}->{rowcol} eq 'row') {
            $sample_cnt = $nrow - 1;
        }
    }
    else {
        my ($nrow, $ncol) = data_check($self->{pic}->{file}, "matrix", "$self->{pic}->{namefix}");
        $sample_cnt = $ncol;
    }

    if($sample_cnt < 1) {
        ERROR("$self->{type}", "no data found, check the data. ");
    }
    else {
        my ($type, $tag, $num) = ("dot", "default", $sample_cnt);
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

    if(defined $self->{pic}->{group} and -e $self->{pic}->{group}) {
        my @return = data_check($self->{pic}->{group}, "group");
        my $group_cnt = $return[1];

        my ($type, $tag, $num) = ("dot", "default", $group_cnt);
        if($self->{pic}->{colset} ne 'none') {
            my @set = split /\_/, $self->{pic}->{colset};
            $type = $set[0];
            $tag  = $set[1];
            $num  = defined $set[2] ? $set[2] : $num;
        }

        my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
        if(scalar(@$colors) >= $group_cnt) {
            $color = join(",", map { "'$_'" } @$colors[0 .. ($group_cnt - 1)]);
        }
        else {
            $color = join(",", map { "'$_'" } @$colors);
            $color = "colorRampPalette(c($color))($group_cnt)";
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
            $color = join(",", map {"'$_'"} @colors);
        }
        else {
            $color = "'$self->{pic}->{collist}'";
        }
    }

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_dot_theme_default();

    my $rscript = <<RSC;
$header
library(gmodels)
#source("/Bio/User/liuyubin/pipeline/GeneralPlot/src/Rlib/all.theme.r")
#mytheme <- dot_theme_default()
$theme

readTable <- function (file, rowcol = "col") {
    data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, 
    row.names = 1, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

file      <- "$self->{pic}->{file}"
group     <- "$self->{pic}->{group}"
arrowfile <- "$self->{pic}->{arrowfile}"
outdir    <- "$self->{pic}->{outdir}"
outname   <- "$self->{pic}->{outname}"

namefix <- "$self->{pic}->{namefix}"
rowcol  <- "$self->{pic}->{rowcol}"

isskipstat  <- as.logical("$self->{pic}->{skipstat}")
pc1         <- as.numeric("$self->{pic}->{pc1}")
pc2         <- as.numeric("$self->{pic}->{pc2}")

title     <- "$self->{pic}->{title}"
xlab      <- "$self->{pic}->{xlab}"
ylab      <- "$self->{pic}->{ylab}"

islegend  <- as.logical("$self->{pic}->{legend}")
islabel   <- as.logical("$self->{pic}->{label}")
isarrow   <- as.logical("$self->{pic}->{arrow}")
isellipse <- as.logical("$self->{pic}->{ellipse}")
isreport  <- as.logical("$self->{pic}->{report}")

dotsize <- "$self->{pic}->{dotsize}"
width   <- $self->{pic}->{width}
height  <- $self->{pic}->{height}
colors  <- c($color)

if(isskipstat == T) {
    pc <- readTable(file, rowcol)
    pc.name <- colnames(pc)
    pc.name.new <- gsub("\\\\(.*?\\\\)", "", pc.name)
    colnames(pc) <- pc.name.new
} else {
    data <- readTable(file, rowcol)
    if(namefix != 'none') {
        getcol <- grep(namefix, colnames(data))
        if(length(getcol) > 0) {
            data <- data[, getcol]
            colnames(data) <- gsub(namefix, "", colnames(data))
        }
    }
    data.mat <- as.matrix(data)
    data.mat.t <- t(data.mat)
    data.pca <- fast.prcomp(data.mat.t, retx = T, scale = F, center = T)
    data.pca.summ <- summary(data.pca)
    data.pca.summ.imp <- data.pca.summ\$importance
    pc <- as.data.frame(data.pca\$x)
    pc.name <- colnames(pc)
}

if(group != "none") {
    classify <- read.table(group, header = F, sep = '\\t', check.names = F, quote = "")
    colnames(classify) <- c("names", "group")
    rownames(classify) <- classify[, 1]
    classify\$group <- factor(classify\$group, levels = unique(classify\$group))
    classify <- classify[rownames(pc), ]
    pc\$group  <- classify\$group
    pc\$names  <- classify\$names
    sample.num <- length(unique(pc\$names))
    group.num  <- length(unique(pc\$group))
    isgroup    <- T
} else {
    pc\$group  <- factor(rownames(pc), levels = rownames(pc))
    pc\$names  <- rownames(pc)
    sample.num <- length(unique(pc\$names))
    group.num  <- length(unique(pc\$group))
    isgroup    <- F
}

if(arrowfile != "none" && isskipstat == T) {
    arrow.r <- readTable(arrowfile, rowcol)
    pc1_id <- colnames(pc)[pc1]
    pc2_id <- colnames(pc)[pc2]
    arrow.d <- data.frame(arrow.r, type = row.names(arrow.r), V_1 = arrow.r[, pc1_id], 
        V_2 = arrow.r[, pc2_id])
    isarrow <- T
}
if(isarrow == T) {
    if(isskipstat != T) {
        pc1_id <- colnames(pc)[pc1]
        pc2_id <- colnames(pc)[pc2]
        pca.rotation <- data.frame(data.pca\$rotation, type = row.names(data.pca\$rotation))
        arrow.d <- transform(pca.rotation, V_1 = get(pc1_id), V_2 = get(pc2_id))
    }
}

if(isskipstat != T) {
    pc.data <- data.frame("id" = rownames(pc), pc)
    pc.pro <- c(as.numeric(sprintf("%.3f", data.pca.summ.imp[2,])) * 100)
    colnames(pc.data)[1:(ncol(pc.data) - 2)] <- c("id", paste(pc.name, "(", pc.pro, "%)", sep = ""))
    if(isreport == T) {
        write.table(pc.data, file = paste0(outdir, "/", outname, ".PC_data.xls"), 
            sep = "\\t", quote = F, row.names = F)
    }
}

if(length(colors) >= group.num) {
    colors <- colors[1:group.num]
} else {
    colors <- colorRampPalette(colors)(group.num)
}

psize = 4
tsize = 4
tvjust = -1.65
if(sample.num <= 20) {
    psize = 6
    tsize = 4
} else if(sample.num <= 40) {
    psize = 5
    tsize = 3.5
} else if(sample.num <= 80) {
    psize = 3
    tsize = 3
} else {
    psize = 2
    tsize = 2
}
if(dotsize != 'none') {
    psize = as.numeric(dotsize)
}

if(isskipstat != T) {
    pro1 <- as.numeric(sprintf("%.3f", data.pca.summ.imp[2,pc1])) * 100
    pro2 <- as.numeric(sprintf("%.3f", data.pca.summ.imp[2,pc2])) * 100
    if(xlab == "none") {
        xlab = paste("PC", pc1, "(", pro1, "%)", sep = "")
    }
    if(ylab == "none") {
        ylab = paste("PC", pc2, "(", pro2, "%)", sep = "")
    }
} else {
    if(xlab == "none") {
        xlab = pc.name[pc1]
    }
    if(ylab == "none") {
        ylab = pc.name[pc2]
    }
}

xmin = min(pc[, pc1]) * 1.2
xmax = max(pc[, pc1]) * 1.2
ymin = min(pc[, pc2]) * 1.2
ymax = max(pc[, pc2]) * 1.2

pc.d <- pc
colnames(pc.d)[pc1] <- "V_1"
colnames(pc.d)[pc2] <- "V_2"
p <- ggplot(pc.d, aes(V_1, V_2))

if(isellipse == T && isgroup == T) {
    dot.sum <- pc.d %>% group_by(group) %>% summarize(count = n())
    if(min(dot.sum\$count) >= 4) {
        p <- p + stat_ellipse(aes(fill = group), type = 'norm', geom  = "polygon",
            alpha = 0.2, color = NA) +
        scale_fill_manual(values = colors)
    }
}

if(isgroup == T) {
    if(group.num > 6) {
        p <- p + geom_point(size = psize, aes(color = group))
    } else {
        point.arr <- c(16, 17, 15, 18, 13, 7)
        point.arr <- point.arr[1:group.num]
        p <- p + geom_point(size = psize, aes(color = group, shape = group)) + 
        scale_shape_manual(values = point.arr)
    }
} else {
    p <- p + geom_point(size = psize, aes(color = group))
}

if(islabel == T) {
    p <- p + geom_text(aes(label = names), size = tsize, vjust = tvjust)
}

if(isarrow == T) {
    arrowshift <- function(x, y, size = 1) {
        dist0 <- sqrt(x^2 + y^2)
        dist  <- dist0 + size
        coeff <- dist / dist0
        coeff
    }
    arrow.max.len <- max(sqrt(arrow.d\$V_1^2 + arrow.d\$V_2^2))
    arrow.d.sh <- arrow.d %>% mutate(coeff = arrowshift(V_1, V_2, size = arrow.max.len / 5), 
        V_1 = V_1 * coeff, V_2 = V_2 * coeff)

    p <- p + geom_segment(data = arrow.d, aes(x = 0, y = 0, xend = V_1, yend = V_2), 
        arrow = arrow(length = unit(0.2, "cm")), size = 0.6, alpha = 0.7, color = 'red')
    p <- p + geom_text(data = arrow.d.sh, aes(V_1, V_2, label = type), size = 4, 
        nudge_x = 0, nudge_y = 0, alpha = 0.7, color = 'red' )
}

p <- p + scale_color_manual(values = colors) + 
scale_x_continuous(breaks = pretty_breaks(5), limits = c(xmin, xmax)) +
scale_y_continuous(breaks = pretty_breaks(5), limits = c(ymin, ymax)) + 
geom_hline(yintercept = 0, linetype = 4, color="grey") +
geom_vline(xintercept = 0, linetype = 4, color="grey") + 
labs(x = xlab, y = ylab, title = title) + 
mytheme + 
guides(color = guide_legend(override.aes = list(size = 4)))

if(islegend == F) {
    p <- p + theme(legend.position = 'none')
} else if(islegend == T) {
    width <- width + ceiling(group.num / 30) * 1.2
}

ggsave(paste0(outdir, "/", outname, ".pdf"), p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "pca";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  general the volcano 
#-------------------------------------------------------------------------------
sub plot_volcano {

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

    $self->{pic}->{header}  = $opts{'-header'}  ? $opts{'-header'} : "T";
    $self->{pic}->{rowcol}  = $opts{'-rowcol'}  ? $opts{'-rowcol'}  : "col";    # row
    $self->{pic}->{namefix} = $opts{'-namefix'} ? $opts{'-namefix'} : "none";   # _fpkm ...

    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";   # 'red,blue'
    
    $self->{pic}->{fccolname} = $opts{'-fccolname'} ? $opts{'-fccolname'} : "log2(fc)";
    $self->{pic}->{pqcolname} = $opts{'-pqcolname'} ? $opts{'-pqcolname'} : "FDR";
    $self->{pic}->{fccut}     = $opts{'-fccut'}     ? $opts{'-fccut'}     : 1;
    $self->{pic}->{pqcut}     = $opts{'-pqcut'}     ? $opts{'-pqcut'}     : 0.05;
    $self->{pic}->{fcformat}  = $opts{'-fcformat'}  ? $opts{'-fcformat'}  : 'none';    # none
    $self->{pic}->{pqformat}  = $opts{'-pqformat'}  ? $opts{'-pqformat'}  : '-log10';  # none

    $self->{pic}->{title}   = defined $opts{'-title'}  ? $opts{'-title'}  : "Volcano Plot";
    $self->{pic}->{xlab}    = defined $opts{'-xlab'}   ? $opts{'-xlab'}   : "log2(fc)";
    $self->{pic}->{ylab}    = defined $opts{'-ylab'}   ? $opts{'-ylab'}   : "-log10(FDR)";
    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 7;

    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "T";

    # color settings
    my $color = "";
    my ($type, $tag, $num) = ("dot", "volcano", "0");
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
    my $theme  = ggplot2_dot_theme_default();

    my $rscript = <<RSC;
$header
$theme

readTable <- function (file, rowcol = "col", row.names = 1) {
    data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, 
    row.names = row.names, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

RSC

    $rscript .= <<RSC;
file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

log2fc_col  <- "$self->{pic}->{fccolname}"
pqvalue_col <- "$self->{pic}->{pqcolname}"
fccut <- $self->{pic}->{fccut}
pqcut <- $self->{pic}->{pqcut}
fcfmt <- "$self->{pic}->{fcformat}"
pqfmt <- "$self->{pic}->{pqformat}"

width <- $self->{pic}->{width}
height <- $self->{pic}->{height}
colors <- c($color)

data  <- readTable(file, row.names = NULL)
fccol <- grep(log2fc_col, colnames(data), fixed = T)
pqcol <- grep(pqvalue_col, colnames(data), fixed = T)
data  <- data[, c(fccol, pqcol)]

colname <- colnames(data)
colnames(data) <- c("fc", "pqvalue")
if(fcfmt == "log2") {
    data[, 1] = log2(data[, 1])
}

## color setting
data\$color = "nosig"
data[ data[,1] > fccut  & data[,2] < pqcut , 3 ] = "up"
data[ data[,1] < -fccut & data[,2] < pqcut , 3 ] = "down"
data\$color <- factor(data\$color, levels = c("up", "nosig", "down"))
colors <- colors[sort(unique(data\$color))]

data[data[,2] == 0, 2] = min(data[data[,2] > 0, 2]) / 10
if(pqfmt == "-log10") {
    data[, 2] = -1 * log10(data[, 2])
}

ymax = max(data[,2]) * 1.1        # quantile(data[,2], probs = c(0.9))
xmax = max(abs(data[,1])) * 1.1   # quantile(abs(data[,1]), probs = c(0.9))

title = "$self->{pic}->{title}"
xlab = colname[1]
ylab = paste0("-log10(", colname[2], ")")

p <- ggplot(data, aes(fc, pqvalue, color = color)) + 
geom_point(size = 0.8, shape = 19) + 
scale_x_continuous(breaks = pretty_breaks(5), limits = c(-xmax, xmax)) +
scale_y_continuous(breaks = pretty_breaks(5), limits = c(0, ymax)) + 
scale_color_manual(values = colors) + 
geom_hline(yintercept = -log10(pqcut), linetype = 2) +
geom_vline(xintercept = c(-fccut, fccut), linetype = 2) +
labs(x = xlab, y = ylab, title = title, color = "") +
mytheme + 
guides(color = guide_legend(override.aes = list(size = 3)))

outfig <- paste0(outdir, "/", outname, ".pdf")
ggsave(outfig, p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "volcano";
    $self->{obj} = 1;

    return $self;

}

#-------------------------------------------------------------------------------
#  general the manhattan 
#-------------------------------------------------------------------------------
sub plot_manhattan {

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
    defined $opts{'-chrlen'}  and $self->{pic}->{chrlen}  = rel2abs($opts{'-chrlen'});
    defined $opts{'-outdir'}  and $self->{pic}->{outdir}  = rel2abs($opts{'-outdir'});
    defined $opts{'-outname'} and $self->{pic}->{outname} = $opts{'-outname'};

    $self->{pic}->{format}  = $opts{'-format'}  ? $opts{'-format'}  : "table";
    $self->{pic}->{datacol} = $opts{'-datacol'} ? $opts{'-datacol'} : "1,2,3,4";
    $self->{pic}->{thres}   = $opts{'-thres'}   ? $opts{'-thres'}   : "0.01,0.05";
    $self->{pic}->{scale}   = $opts{'-scale'}   ? $opts{'-scale'}   : "-log10";   # -log2 | -log10
    $self->{pic}->{space}   = $opts{'-space'}   ? $opts{'-space'}   : "T";
    $self->{pic}->{dotsize} = $opts{'-dotsize'} ? $opts{'-dotsize'} : 1.3;

    $self->{pic}->{colset}  = $opts{'-colset'}  ? $opts{'-colset'}  : "none";
    $self->{pic}->{colnpg}  = $opts{'-colnpg'}  ? $opts{'-colnpg'}  : "F";
    $self->{pic}->{collanc} = $opts{'-collanc'} ? $opts{'-collanc'} : "F";
    $self->{pic}->{collist} = $opts{'-collist'} ? $opts{'-collist'} : "none";     # 'red,blue'

    $self->{pic}->{header}  = $opts{'-header'} ? $opts{'-header'} : "T";
    $self->{pic}->{title}   = defined $opts{'-title'}  ? $opts{'-title'}  : "Manhattan Plot";
    $self->{pic}->{xlab}    = defined $opts{'-xlab'}   ? $opts{'-xlab'}   : "Chromosomes";
    $self->{pic}->{ylab}    = defined $opts{'-ylab'}   ? $opts{'-ylab'}   : "-log10(Pvalue)";
    $self->{pic}->{legend}  = $opts{'-legend'} ? $opts{'-legend'} : "F";

    $self->{pic}->{width}   = $opts{'-width'}  ? $opts{'-width'}  : 8;
    $self->{pic}->{height}  = $opts{'-height'} ? $opts{'-height'} : 6;

    # color settings
    my $label_cnt = 0;
    #-------------------------------------
    # we assume the data format is like: 
    # SNP   CHR POS P
    # rs1   1   1   0.9148060
    # rs2   1   2   0.9370754
    # rs3   1   3   0.2861395
    my @return = data_check($self->{pic}->{file}, "table2");
    if(!defined $return[1]) {
        ERROR("plot_manhattan", "check the input data . ");
    }
    $label_cnt = $return[1];

    my $color = "";
    my ($type, $tag, $num) = ("dot", "default", $label_cnt);
    if($self->{pic}->{colset} ne 'none') {
        my @set = split /\_/, $self->{pic}->{colset};
        $type = $set[0];
        $tag  = $set[1];
        $num  = defined $set[2] ? $set[2] : $num;
    }

    my $colors = fetch_color(-type => $type, -tag => $tag, -num => $num);
    if(scalar(@$colors) >= $label_cnt) {
        $color = join(",", map {"'$_'"} @$colors[0 .. ($label_cnt - 1)]);
    }
    else {
        $color = join(",", map {"'$_'"} @$colors);
        $color = "colorRampPalette(c($color))($label_cnt)";
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
        my @colors = split /,/, $self->{pic}->{collist};
        $color = join(",", map {"'$_'"} @colors);
    }

    # theme settings
    my $header = ggplot2_head();
    my $theme  = ggplot2_dot_theme_manhattan();

    my $rscript = <<RSC;
$header
$theme
source("$SOFT{ggplot2}")

RSC

    if(defined $self->{pic}->{chrlen} and -e $self->{pic}->{chrlen}) {
        $rscript .= <<RSC;
file    <- "$self->{pic}->{file}"
chrlen  <- "$self->{pic}->{chrlen}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

RSC
    }
    else {
        $rscript .= <<RSC;
file    <- "$self->{pic}->{file}"
outdir  <- "$self->{pic}->{outdir}"
outname <- "$self->{pic}->{outname}"

RSC
    }

    $rscript .= <<RSC;
datacol <- as.numeric(c($self->{pic}->{datacol}))
thres   <- as.numeric(c($self->{pic}->{thres}))
scale   <- "$self->{pic}->{scale}"
space   <- as.logical("$self->{pic}->{space}")
dotsize <- as.numeric("$self->{pic}->{dotsize}")

title   <- "$self->{pic}->{title}"
xlab    <- "$self->{pic}->{xlab}"
ylab    <- "$self->{pic}->{ylab}"
legend  <- as.logical("$self->{pic}->{legend}")

width   <- $self->{pic}->{width}
height  <- $self->{pic}->{height}
colors <- c($color)

if(length(datacol) != 4 || is.numeric(datacol) != T) {
    cat("Error: you should set 4 number of columns for manhattan data: [SNP CHR POS P]. \\n")
    q()
}
if(is.na(thres) == T || is.numeric(thres) != T) {
    cat("Warn: the threshold [", thres, "] is irregular and would not be applied. \\n")
    thres_tag = F
} else {
    thres_tag = T
}

data <- read.table(file, head = $self->{pic}->{header}, sep = "\\t", check.names = F, quote = "")
data <- data[, datacol]
colnames(data) <- c("SNP", "CHR", "POS", "P")

if(scale == "log10") {
    data\$logP <- log10(data\$P)
} else if(scale == "log2") {
    data\$logP <- log2(data\$P)
} else if(scale == "-log10") {
    data\$logP <- -log10(data\$P)
} else if(scale == "-log2") {
    data\$logP <- -log2(data\$P)
} else {
    data\$logP <- data\$P
}

RSC

    if(defined $self->{pic}->{chrlen} and -e $self->{pic}->{chrlen}) {
        $rscript .= <<RSC;
chr_len <- read.table(chrlen, header = F, sep = '\\t', check.names = F, quote = "")
chr_len <- as.tbl(chr_len[, 1:2])
colnames(chr_len) <- c("CHR", "chr_len")

RSC
    }
    else {
        $rscript .= <<RSC;
chr_len <- data %>% group_by(CHR) %>% summarise(chr_len = max(POS))

RSC
    }

    $rscript .= <<RSC;
chr_pos <- chr_len %>% mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>% select(-chr_len)
if(space == T) {
    space = round(mean(chr_len\$chr_len) / 4)
    chr_pos <- chr_len %>% mutate(total = cumsum(as.numeric(chr_len)) - chr_len, space = space) %>% 
    mutate(cumspace = cumsum(space) - space, total_space = total + cumspace) %>% 
    select(CHR, total_space)
    colnames(chr_pos)[2] = "total"
}

snp_pos <- chr_pos %>% left_join(data, ., by="CHR") %>% arrange(CHR, POS) %>% 
    mutate( POScum = POS + total)

# filter NA data  
snp_pos <- snp_pos %>% filter(is.na(POScum) != T)
x_axis  <- snp_pos %>% group_by(CHR) %>% summarize(center=(max(POScum) + min(POScum)) / 2)

chr_cnt <- length(unique(x_axis\$CHR))
if(chr_cnt != length(colors)) {
    colors <- colorRampPalette(colors)(chr_cnt)
}
#colors <- rep(c("grey", "skyblue"), int(chr_cnt / 2 + 1))

if(chr_cnt <= 20) {
    width = 8
} else if(chr_cnt <= 30) {
    width = 10
} else if(chr_cnt <= 100) {
    width = 10 + 8 * (chr_cnt - 30) / 70
} else {
    width = 20
}
axis.size <- AxisTextSize(x_axis\$CHR, width = width)

# recheck for min distance of chrs
#chr_sumlen <- sum(chr_len\$chr_len)
chr_sumlen <- 2 * x_axis[chr_cnt, 2] - x_axis[chr_cnt - 1, 2]
min_dist <- min(abs(x_axis\$center[2:chr_cnt] - x_axis\$center[1:(chr_cnt-1)]))

x_angle_tag = F
min_ratio <- min_dist / chr_sumlen
if(min_ratio < 0.025) {
    tmp_ratio <- 0.02
    axis.size <- width * tmp_ratio * 72 * 0.7
}
if(min_ratio < 0.025) {
    x_angle_tag = T
}

ymin <- 0
ymax <- max(snp_pos\$logP) * 1.3
colors.thres <- c("red", "green", "blue")
if(thres_tag == T) {
    thres <- sort(thres)
    if(scale == "-log10") {
        thres <- -log10(thres / nrow(snp_pos))
    } else if(scale == "-log2") {
        thres <- -log2(thres / nrow(snp_pos))
    }
    if(length(colors.thres) != length(thres)) {
        colors.thres <- colorRampPalette(colors.thres)(length(thres))
    }
    ymax <- max(thres) * 1.2
}

p <- ggplot(snp_pos, aes(x = POScum, y = logP)) + 
geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = dotsize) + 
scale_color_manual(values = colors) + 
scale_x_continuous(label = x_axis\$CHR, breaks = x_axis\$center, expand = c(0.01, 0.01)) + 
scale_y_continuous(expand = c(0.01, 0.01), limits = c(ymin, ymax)) + 
labs(x = xlab, y = ylab, title = title) +
mytheme + 
theme(axis.text = element_text(size = axis.size)) 

if(x_angle_tag == T) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}
if(thres_tag == T) {
    p <- p + geom_hline(yintercept = thres, size = 0.8, color = colors.thres, 
        linetype = rep("twodash", length(thres)), alpha = 0.8) 
}
if(legend == F) {
    p <- p + theme(legend.position = 'none')
}

ggsave(paste0(outdir, "/", outname, ".pdf"), p, width = width, height = height)
ggsave(paste0(outdir, "/", outname, ".png"), p, width = width, height = height)

RSC

    $self->{rscript} = $rscript;
    $self->{method} = "manhattan";
    $self->{obj} = 1;

    return $self;

}

return 1;

