package GeneralPlot::Help;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.05.22 15:07:55          |
#--------------------------------------------------#
use strict;
use warnings;
no strict 'refs';
no warnings 'qw';
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Config::General;
use parent 'Exporter';
our @EXPORT = qw(type_check usage sub_usage);
our @EXPORT_OK = qw();

use lib "$Bin";
use lib "$Bin/..";

use GeneralPlot::Debug;

sub type_check {

    my $type = shift;

    my @existed_type = (qw/venn scatter pca volcano manhattan line 
    bar errorbar histgram density pie twopie box violin splitviolin 
    correlation correlation_old cluster scatter3d polarbar sankey 
    scatter_bubble lollipop cumulative_line/);

    my %check = ();
    map { $check{$_} = 1 } @existed_type;
    
    if(defined $check{$type}) {
        return 0;
    }
    else {
        return 1;
    }

}

sub usage {

    print <<USAGE;

Usage: 

    perl $0 <GRAPH> [options]

    the valid graphs are:

    venn
    scatter
    scatter3d
    scatter_bubble
    lollipop
    pca
    volcano
    manhattan
    line
    cumulative_line
    bar
    errorbar
    polarbar
    sankey
    histgram
    density
    pie
    twopie
    box
    violin
    splitviolin
    correlation
    correlation_old
    cluster

USAGE

}

sub sub_usage {
    
    my $type = shift;

    my $sub = "_$type";
    my $msg = $sub->();

    print $msg;

    return 0;
}

sub _scatter3d {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Dot->plot(%opts)

    -method         => 'scatter3d'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -colnpg         => F
    -collanc        => F
    -collist        => 'none'

    -x              => 2
    -y              => 3
    -z              => 4
    -samplecol      => 1
    -groupcol       => none

    -header         => T
    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -zlab           => ''
    -pcex           => 1
    -tcex           => 0.8
    -vjust          => 0.5

    -grid           => F
    -label          => F

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
   
    $usage;

}

sub _cumulative_line {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Line->plot(%opts)

    -method         => 'cumulative_line'

    necessary:
    -file | -f
    -outprefix | -op
    -outdir
    -outname

    optional:
    -group          => 'none'
    -scale          => 'none'   # none | log2 | log10 
    -filter0        => F | T
    -namefix        => 'none'   # _fpkm ...
    -color          => 'none'   # 'red,blue'
    -title          => 'none'
    -x_lab          => 'none'
    -y_lab          => 'none'
    -size           => 1        # line size
    -alpha          => 0.9      # line transparency
    -width          => 7
    -height         => 7
    -out_format     => 'pic_data'   # pic_data | pic | data

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
   
    $usage;

}

sub _lollipop {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Dot->plot(%opts)

    -method         => 'lollipop'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -group          => 'none'
    -scale          => 'none'   # none | log2 | log10 | '-log2' | '-log10'
    -dot_size       => 4
    -line_size      => 1
    -dot_color      => "none'
    -line_color     => "grey50" # none means no line
    -transparency   => 1
    -is_label       => F
    -is_baseline    => F
    -is_flip        => F        # coord flip
    -is_sort        => F        # sort the data by group and value
    -title          => 'none'
    -x_lab          => 'none'
    -y_lab          => 'none'
    -lgd_lab        => 'none'

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
   
    $usage;

}


sub _scatter_bubble {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Dot->plot(%opts)

    -method         => 'scatter_bubble'

    necessary:
    -file1          # exp_file
    -file2          # pvalue_file
    -outprefix
    -outdir
    -outname

    optional:
    -rowname_file   => 'none'
    -colname_file   => 'none'
    -scale1         => 'none'   # none | log2 | log10 | '-log2' | '-log10'
    -scale2         => 'none'   # none | '-log2' | '-log10' | log2 | log10
    -dot_minsize    => 1
    -dot_maxsize    => 6
    -transparency   => 0.9
    -collist        => 'none'   # 'red,blue'
    -title          => 'none'
    -xlab           => 'none'
    -ylab           => 'none'
    -lgd1           => 'none'   # default: none => value.1
    -lgd2           => 'none'   # default: none => value.2

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
   
    $usage;

}

sub _venn {

    my $usage = <<USAGE;

    GeneralPlot::Image::Venn->plot(%opts)

    -method         => 'venn'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -width          => 600
    -height         => 600
    -margin         => 30
    -inwidth        => 1

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _scatter {

    my $usage = <<USAGE;

    GeneralPlot::Image::Dot->plot(%opts)

    -method         => 'dot'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -colnpg         => F
    -collanc        => F
    -collist        => 'none'

    -x              => 1
    -y              => 2
    -group          => none

    -header         => T
    -title          => "Scatter Plot"
    -xlab           => ''
    -ylab           => ''
    -legend         => T
    -width          => 8
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
   
    $usage;

}

sub _pca {

    my $usage = <<USAGE;

    GeneralPlot::Image::Dot->plot(%opts)

    -method         => 'pca'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -group
    -arrowfile

    -skipstat       => F

    -pc1            => 1
    -pc2            => 2
    
    -colnpg         => F
    -collanc        => F
    -collist        => 'none'

    -arrow          => F
    -ellipse        => F
    -report         => T

    -header         => T
    -rowcol         => 'col'    [ col | row ]
    -namefix        => 'none'   [ _fpkm ... ]

    -title          => "PCA"
    -label          => F
    -legend         => T
    -dotsize        => 'none'
    -width          => 8
    -height         => 8

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _volcano {

    my $usage = <<USAGE;

    GeneralPlot::Image::Dot->plot(%opts)

    -method         => 'volcano'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -rowcol         => 'col'    [ col | row ]
    -namefix        => 'none'   [ _fpkm ... ]

    -fccolname      => 'log2(fc)'
    -pqcolname      => 'FDR'
    -fccut          => 1
    -pqcut          => 0.05
    -fcformat       => 'none'     [ none | log2 ]
    -pqformat       => '-log10'   [ none | -log10 ]

    -collist        => 'none'

    -title          => "Volcano Plot"
    -xlab           => 'log2(fc)'
    -ylab           => '-log10(FDR)'
    -width          => 8
    -height         => 7
    -legend         => "T"

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _manhattan {

    my $usage = <<USAGE;

    GeneralPlot::Image::Dot->plot(%opts)

    -method         => 'manhattan'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -chrlen

    -header         => T
    -datacol        => '1,2,3,4'
    -thres          => '0.01,0.05',
    -scale          => '-log10',
    -space          => 'T',

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -title          => "Manhattan Plot"
    -xlab           => 'Chromosomes'
    -ylab           => '-log10(Pvalue)'
    -legend         => "F"

    -width          => 8
    -height         => 6

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _line {

    my $usage = <<USAGE;

    GeneralPlot::Image::Line->plot(%opts)

    -method         => 'line'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -legend         => "T"

    -width          => 7
    -height         => 6

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _stauration {

    my $usage = <<USAGE;

    GeneralPlot::Image::Line->plot(%opts)

    -method         => ''

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _bar {

    my $usage = <<USAGE;

    GeneralPlot::Image::Bar->plot(%opts)

    -method         => 'bar'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -barfill        => 'F'
    -bardodge       => 'F'
    -mulcolor       => 'F'
    -bracket        => 'T'

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',
    -coordflip      => 'F',

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -label          => 'F'
    -legend         => 'T'

    -width          => 8
    -height         => 6

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _errorbar {

    my $usage = <<USAGE;

    GeneralPlot::Image::Bar->plot(%opts)

    -method         => 'errorbar'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -mulcolor       => 'F'

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -label          => 'F'
    -legend         => 'T'

    -width          => 8
    -height         => 6

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _polarbar {

    my $usage = <<USAGE;

    GeneralPlot::Image::Bar->plot(%opts)

    -method         => 'polarbar'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -collist        => 'none'
    -colset         => 'rainbow'    # rainbow | viridis | greyred
    -title          => 'none'
    -format         => 'none'       # log2 | log10 
    -header         => T
    -show_border    => F
    -show_label     => T
    -show_count     => T
    -direction      => 'clockwise'  # clockwise | anticlockwise
    -display        => 'type1'      # type1 | type2
    -width          => 8
    -height         => 8

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _sankey {

    my $usage = <<USAGE;

    GeneralPlot::Image::Bar->plot(%opts)

    -method         => 'sankey'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -collist        => 'none'
    -colset         => 'set1'       # set1 set2 set3 set4
    -transparency   => 0.7          # 0 ~ 1
    -band_width     => 0.15         # 0 ~ 1
    -band_knotpos   => 0.4          # 0 ~ 1
    -title          => 'none'
    -xlab           => 'none'
    -ylab           => 'none'
    -title_size     => 16
    -axis_size      => 12
    -label_size     => 4
    -label_color    => "black"
    -stratum_color  => "black"
    -alluvium_color => NA
    -fill_column    => "default"    # can be default(the first column), or number (1, 2, 3 ...)

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _histgram {

    my $usage = <<USAGE;

    GeneralPlot::Image::Hist->plot(%opts)

    -method         => 'hist'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -format         => 'table'      [ table | matrix ]
    -colname        => 'NULL'

    -facet          => 'F'
    -identity       => 'F'
    -density        => 'F'

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'

    -width          => 8
    -height         => 8
    -bins           => 30

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE
    
    $usage;

}

sub _density {

    my $usage = <<USAGE;

    GeneralPlot::Image::Hist->plot(%opts)

    -method         => 'density'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -format         => 'table'      [ table | matrix ]
    -colname        => 'NULL'

    -namefix        => 'none'       # _fpkm ...
    -scale          => "none"       # log2 | log10
    -filter0        => 'F' 

    -facet          => 'F'
    -inline         => 'F'

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'

    -width          => 7
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _pie {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Pie->plot(%opts)

    -method         => 'pie'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -direction      => -1      # clockwise(1) | anticlockwise(-1)
    -order          => 'none'  # none | ab | ba

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -label          => 'F'
    -legend         => 'T'

    -width          => 7
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _twopie {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Pie->plot(%opts)

    -method         => 'twopie'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -label          => 'F'
    -legend         => 'T'

    -width          => 7
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _box {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Box->plot(%opts)

    -method         => 'box'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -format         => 'matrix'     [ matrix | table ]
    -rowcol         => 'col'        [ col | row ]
    -group          => 'none'

    -namefix        => 'none'       [ _fpkm ... ]
    -filter0        => 'F'
    -scale          => 'none'       [ log2 | log10 ]
    -ymin           => -999999
    -ymax           => 999999

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -report         => 'T'

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'
    -adddot         => 'F'

    -width          => 8
    -height         => 6

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _violin {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Box->plot(%opts)

    -method         => 'violin'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -format         => 'matrix'     [ matrix | table ]
    -rowcol         => 'col'        [ col | row ]
    -group          => 'none'

    -namefix        => 'none'       [ _fpkm ... ]
    -filter0        => 'F'
    -scale          => 'none'       [ log2 | log10 ]
    -ymin           => -999999
    -ymax           => 999999

    -colnpg         => 'F',
    -collanc        => 'F',
    -collist        => 'none',

    -report         => 'T'

    -title          => "Violin Plot"
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'
    -adddot         => 'F'

    -width          => 8
    -height         => 6

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _splitviolin {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Box->plot(%opts)
    Notice: can only be used on 2 attrs data

    -method         => 'splitviolin'

    necessary:
    -file
    or:
    -file1
    -file2
    -name1
    -name2

    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -rowcol         => 'col'        [ col | row ]
    -namefix        => 'none'       [_fpkm ...]
    -filter0        => 'F'
    -scale          => 'none'       [ log2 | log10 ]
    -ymin           => -999999
    -ymax           => 999999

    -collist        => 'none',
    -alpha          => 0.75,
    -bartype        => 'split'      [split | center | none]

    -report         => 'T'

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'
    -adddot         => 'F'

    -width          => 8
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _heatmap {

}

sub _correlation {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Heatmap->plot(%opts)

    -method         => 'ggcor'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -rowcol         => 'col'        [ col | row ]
    -namefix        => 'none'       [ _fpkm ... ]

    -cor            => 'pearson'    [ pearson | spearman | kendall ]
    -cortest        => 'F'
    -pvalue         => 0.05

    -postype        => 'full'       [ upper | lower ]

    -collist        => 'none',

    -title          => "Sample Correlation"
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'
    -label          => 'F'

    -width          => 7
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _correlation_old {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Heatmap->plot(%opts)

    -method         => 'corr'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -rowcol         => 'col'        [ col | row ]
    -namefix        => 'none'       [ _fpkm ... ]

    -cor            => 'pearson'    [ pearson | spearman | kendall ]

    -collist        => 'none',

    -title          => "Sample Correlation"
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'
    -label          => 'F'

    -width          => 7
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

sub _cluster {

    my $usage = <<USAGE;
    
    GeneralPlot::Image::Heatmap->plot(%opts)

    -method         => 'cluster'

    necessary:
    -file
    -outprefix
    -outdir
    -outname

    optional:
    -header         => T
    -rowcol         => 'col'        [ col | row ]
    -namefix        => 'none'       [ _fpkm ... ]

    -clustrow       => T,
    -clustcol       => T,

    -collist        => 'none',

    -title          => ""
    -xlab           => ''
    -ylab           => ''
    -legend         => 'T'
    -label          => 'F'

    -width          => 7
    -height         => 7

    -save           => T
    -run            => T
    -clean          => T
    -log            => F

USAGE

    $usage;

}

return 1;

