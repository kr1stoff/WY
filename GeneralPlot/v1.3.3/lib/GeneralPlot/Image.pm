package GeneralPlot::Image;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.06.04 05:15:21          |
#--------------------------------------------------#
use strict;
use warnings;
no strict 'refs';
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Config::General;
use parent 'Exporter';
our @EXPORT = qw();
our @EXPORT_OK = qw();

use lib "$RealBin";
use lib "$RealBin/..";

use GeneralPlot::Envs;
use GeneralPlot::Debug;

use GeneralPlot::Image::Venn;
use GeneralPlot::Image::Line;
use GeneralPlot::Image::Pie;
use GeneralPlot::Image::Bar;
use GeneralPlot::Image::Hist;
use GeneralPlot::Image::Box;
use GeneralPlot::Image::Dot;
use GeneralPlot::Image::Heatmap;

#-------------------------------------------------------------------------------
#  create a object 
#-------------------------------------------------------------------------------
sub new {

    my ($class, %opts) = @_;

    my $self = {
        'type'   => 'none',
        'method' => 'none',
        'obj'    => 0,
        'save'   => 0,
        'run'    => 0,
    };

    init($self, %opts);

    bless $self, $class;

    return $self;

}

sub init {

    my ($self, %opts) = @_;

    for my $opt (keys %opts) {
        $self->{$opt} = $opts{$opt};
    }

    return $self;

}

#-------------------------------------------------------------------------------
# plot
#-------------------------------------------------------------------------------
sub plot {

    my ($self, %opts) = @_;

    my %types = (

        'venn'              => 'venn',
        
        'scatter'           => 'dot',
        'scatter3d'         => 'dot',
        'pca'               => 'dot',
        'volcano'           => 'dot',
        'manhattan'         => 'dot',
        'scatter_bubble'    => 'dot',
        'lollipop'          => 'dot',
        
        'line'              => 'line',
       #'stauration'        => 'line',
        'cumulative_line'   => 'line',

        'bar'               => 'bar',
        'errorbar'          => 'bar',
        'polarbar'          => 'bar',
        'sankey'            => 'bar',

        'histgram'          => 'hist',
        'density'           => 'hist',

        'pie'               => 'pie',
        'twopie'            => 'pie',

        'box'               => 'box',
        'violin'            => 'box',
        'splitviolin'       => 'box',

       #'heatmap'           => 'heatmap',
        'correlation'       => 'heatmap',
        'correlation_old'   => 'heatmap',
        'cluster'           => 'heatmap',

    );

    my %methods = (

        'venn'              => 'venn',     
        
        'scatter'           => 'dot',
        'scatter3d'         => 'scatter3d',
        'pca'               => 'pca',
        'volcano'           => 'volcano',
        'manhattan'         => 'manhattan',
        'scatter_bubble'    => 'scatter_bubble',
        'lollipop'          => 'lollipop',

        'line'              => 'line',
        'cumulative_line'   => 'cumulative_line',
       #'stauration'        => 'stauration',

        'bar'               => 'bar',
        'errorbar'          => 'errorbar',
        'polarbar'          => 'polarbar',
        'sankey'            => 'sankey',

        'histgram'          => 'hist',
        'density'           => 'density',

        'pie'               => 'pie',
        'twopie'            => 'twopie',

        'box'               => 'box',
        'violin'            => 'violin',
        'splitviolin'       => 'splitviolin',

        #'heatmap'          => 'heatmap',
        'correlation'       => 'ggcor',
        'correlation_old'   => 'corr',
        'cluster'           => 'cluster',

    );

    $opts{'-save'}  ||= 'T';
    $opts{'-run'}   ||= 'T';
    $opts{'-clean'} ||= 'T';
    $opts{'-log'}   ||= 'F';

    # a abbr for input/iutput
    defined $opts{'-f'}  and !defined $opts{'-file'} and $opts{'-file'} = $opts{'-f'};
    defined $opts{'-op'} and !defined $opts{'-outprefix'} and $opts{'-outprefix'} = $opts{'-op'};

    if(defined $opts{'-graph'}) {
        $self->{type} = defined $types{$opts{'-graph'}} ? 
        $types{$opts{'-graph'}} : $self->{type};
        
        $self->{method} = defined $methods{$opts{'-graph'}} ? 
        $methods{$opts{'-graph'}} : $self->{method};
    }

    if($self->{type} eq 'none' or $self->{method} eq 'none') {
        ERROR('plot', "<-graph> [$opts{'-graph'}] should be set rightly. ");
    }

    if(defined $opts{'-log'} and $opts{'-log'} =~ /^T$|^True$|1/i) {
        unlink("$opts{'-outprefix'}.run.log") if(-e "$opts{'-outprefix'}.run.log");
        open STDOUT, ">> $opts{'-outprefix'}.run.log" or die "[Error]: cannot write to log \n";
        open STDERR, ">> $opts{'-outprefix'}.run.log" or die "[Error]: cannot write to log \n";
    }

    my $plot = launch($self->{type}, $self->{method}, %opts);
    if(defined $opts{'-save'} and $opts{'-save'} =~ /^T$|^True$|1/i) {
        $plot->save(%opts);
    }
    if(defined $opts{'-run'} and $opts{'-run'} =~ /^T$|^True$|1/i) {
        $plot->run(%opts);
    }

    return $plot;

}

sub launch {

    my ($type, $method, %opts) = @_;

    $type = lc($type);
    $method = lc($method);
    my $module = "GeneralPlot::Image::" . ucfirst($type);

    my $self = $module->new();
    $self->plot(-method => $method, %opts);

    return $self;
}

return 1;

