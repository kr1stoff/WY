package SBV::Coordinate::GENE;
#-----------------------------------------------+
#    [APM] This moudle was created by amp.pl    |
#    [APM] Created time: 2018-11-17 11:07:56    |
#-----------------------------------------------+
=pod

=head2 v1.0

Date: 2018-11-17 11:07:56

=head1 Name

SBV::Coordinate::GENE

=head1 Synopsis

This module is not meant to be used directly

=head1 Feedback

Author: Peng Ai
Email:  aipeng0520@163.com

=head1 Version

Version history

=cut


use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT    = qw( );
our @EXPORT_OK = qw( );

use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/../";
use lib "$FindBin::RealBin/../lib";

use SBV::DEBUG;

#===  FUNCTION  ================================================================
#         NAME: new
#      PURPOSE: init the GENE coordinate
#   PARAMETERS: -exons, [-orient, -start, -end, -intron_fix, -intron_scale]
#      RETURNS: a GENE coordinate object
#  DESCRIPTION: ????
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub new {
    my $class = shift ;

    my %opts = @_;
    my $gene = {};

    ERROR("no_exons") unless $opts{'-exons'};
    
    # reorder the input exons
    my @exons = sort { $_->[0] > $_->[1] ? [$_->[1],$_->[0]] : $_ } sort { $a->[0] <=> $b->[0] } @{$opts{'-exons'}};
    
    # set the attributes of gene coordinate
    $gene{exons} = \@exons;

    if ($opts{'-orient'} eq "+") {
        $gene{start} = $exons[0][0];
        $gene{end}   = $exons[-1][1];
    } else {
        $gene{start} = $exons[-1][1];
        $gene{end}   = $exons[0][0];
    }
    
    

    bless $gene , $class;
    return $gene;
}

sub parent {
    my $self = shift;
    $self->{parent} = shift if (@_);
    return $self->{parent};
}


#-------------------------------------------------------------------------------
#  set the size of Gene axis
#  attrs: x1, x2 ,y1, y2
#-------------------------------------------------------------------------------
sub size {
    my $self = shift
    my %opts = @_;
}
