package GeneralPlot::Image::Venn;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.10 16:12:32          |
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

use GeneralPlot::Debug;
use GeneralPlot::Colors;
use GeneralPlot::Data;
use GeneralPlot::Theme;
use GeneralPlot::Envs;

#-------------------------------------------------------------------------------
#  Implementation by SBV
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#  new a venn object 
#-------------------------------------------------------------------------------
sub new {

    my ($class, %opts) = @_;

    my $self = {
        'type'   => 'venn',
        'method' => 'venn',
        'obj'    => 0,
        'save'   => 0,
        'run'    => 0,
    };
    bless $self, $class;

    return $self;

}

sub plot {

    my ($self, %opts) = @_;

    $opts{'-method'} ||= "venn";

    if(defined $opts{'-method'} and $opts{'-method'} eq 'venn') {
        $self->{method} = $opts{'-method'};
        venn($self, %opts);
    }
    else {
        ERROR("$self->{type}", "check the -method [$opts{'-method'}]. ");
    }

    return $self;

}

#-------------------------------------------------------------------------------
#  general sbv venn conf file
#  options: 
#  -file, -outdir, -outname, -format, -width, -height, -margin, -inwidth
#-------------------------------------------------------------------------------
sub venn {

    my ($self, %opts) = @_;

    if(defined $opts{'-outprefix'}) {
        $opts{'-outdir'}  ||= dirname($opts{'-outprefix'});
        $opts{'-outname'} ||= basename($opts{'-outprefix'});
    }

    !defined $opts{'-file'}     and ERROR("$self->{type}", "<-file> is not defined. ");
    !defined $opts{'-outdir'}   and ERROR("$self->{type}", "<-outdir> is not defined. "); 

    $opts{'-file'} = rel2abs($opts{'-file'});
    $opts{'-outdir'} = rel2abs($opts{'-outdir'});

    $self->{pic}->{file}    = $opts{'-file'};
    $self->{pic}->{outdir}  = $opts{'-outdir'};
    $self->{pic}->{outname} = $opts{'-outname'} ? $opts{'-outname'} : "sbv.venn";

    $self->{pic}->{format}  = $opts{'-format'}  ? $opts{'-format'}  : "list2";
    $self->{pic}->{width}   = $opts{'-width'}   ? $opts{'-width'}   : 600;
    $self->{pic}->{height}  = $opts{'-height'}  ? $opts{'-height'}  : 600;
    $self->{pic}->{margin}  = $opts{'-margin'}  ? $opts{'-margin'}  : 30;
    $self->{pic}->{inwidth} = $opts{'-inwidth'} ? $opts{'-inwidth'} : 1;
    
    ## count the samples, the venn data should be in the format list2. no format check here 
    my $count = 0;
    if(-f $self->{pic}->{file}) {
        ($count) = data_check($self->{pic}->{file}, $self->{pic}->{format});
        $count < 2 and ERROR("$self->{type}", "the number of samples in data is less than 2.");
    }
    else {
        WARN("the file <$self->{pic}->{file}> is not existed. ");
    }

    ## theme settings
    my $theme = "";
    if($count == 4) {
        $self->{pic}->{width}   = $opts{'-width'}   ? $opts{'-width'}   : 600;
        $self->{pic}->{height}  = $opts{'-height'}  ? $opts{'-height'}  : 360;

        $theme = sbv_venn_theme_default(
            '-text_font-size' => 26,
            '-text.label_font-size' => 28,
            '-ellipse_stroke-width' => 0.1, 
            '-ellipse_stroke' => '\#666666',
            '-ellipse_fill-opacity' => 0.5,
        );
    }
    elsif($count == 5) {
        $theme = sbv_venn_theme_default(
            '-path_stroke-width' => 0.1,
            '-path_stroke' => '\#666666',
            '-path_fill-opacity' => 0.5,
        );
    }
    else {
        $theme = sbv_venn_theme_default();
    }

    ## color settings
    my $color = "";
    if($count > 7) {
        my $colors = fetch_color('-type' => 'venn', '-num' => 0);
        my $colors_str = join(",", @$colors);
        my $colors_par = `$SOFT{Rscript} $SOFT{get_color_pal} '$colors_str' $count`;
        $color = join(",", map { "\\$_" } split /,/, $colors_par);
    }
    else {
        my $colors = fetch_color('-type' => 'venn', '-num' => $count);
        $color = join(",", map { "\\$_" } @$colors);
    }

    ## generate the conf file
    my $conf = <<CONF;
# background
width = $self->{pic}->{width}
height = $self->{pic}->{height}
margin = $self->{pic}->{margin}

dir = $self->{pic}->{outdir}
file = $self->{pic}->{outname}

$theme

<venn>

width = $self->{pic}->{inwidth}
file = $self->{pic}->{file}

# the format of input file, default list2
# list2 like :
# Group1	v1,v2,v3
# Group2	v1,v3,v4,v5
#------------------------
# list3 like :
# Group1	Group2
# v1	v1
# v2	v3
# v5	v6
#		v8
format = $self->{pic}->{format}

# the color of the venn fill
# default is CMG_Lee (the styles of the map in wiki)
# rainbow
# none : not fill color
fill = $color

show_label = yes
show_logical_label = no
print_stat_info = yes
stat_file = $self->{pic}->{outdir}/$self->{pic}->{outname}.stat.txt

# set the venn model 
# just for 5 sets and 2 sets
# for 5sets, normal or Branko will create the 5 sets venn figue like Branko (5 ellipses),
#            otherwise will create the 5 sets venn figure like 5 Pear-Shaped figure, default
# for 2 sets, horizontal means two horizontal circles, 
#             vertical means two vertical circles, default
# model = ellipse

</venn>

CONF

    #$self->{rm_file}->{stat} = "$self->{pic}->{outdir}/$self->{pic}->{outname}.stat.txt";

    $self->{conf} = $conf;
    $self->{method} = "venn";
    $self->{obj} = 1;

    return $self;

}

return 1;

