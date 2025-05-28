#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2019-10-11 14:47:00    |
#-----------------------------------------------+
# name: drawColPaletters.pl
# func: 
# version: 1.0

use strict;
use warnings;

use SVG;
use General qw/hex_rgb/;

my $file = shift @ARGV;
my %conf = read_conf($file);
my @pals = fetch_pals(keys %conf);
my $pals_num = scalar @pals;

my $name_w = 120;
my $unit_w = 40;
my $unit_h = 20;
my $spacing = 6;
my $margin = 20;

my $width  = $margin*2 + $name_w + $unit_w*20;
my $height = $margin*2 + $unit_h*$pals_num + $spacing*($pals_num-1);

my $svg = SVG->new(width=>$width,height=>$height);
my $text_style = "font-weight:bold;font-size:12px;dominant-baseline:middle;font:Arial";

my $y = $margin;
for my $i ( 0 .. $#pals ){
    my $x = $margin;
    my $pal_name = $pals[$i];
    
    my $palg = $svg->group(id=>$pal_name);

    # draw name 
    $palg->text(x=>$x+$name_w - $spacing,y=>$y+$unit_h/2,style=>$text_style,'text-anchor'=>"end")->cdata($pal_name);
    
    # draw rect box
    $x += $name_w;
    my $num = (split /-/,$pal_name)[1];
    my @colors = map { $conf{"${pal_name}-$_"} } 1 .. $num;
    
    for (@colors){
        $palg->rect(x=>$x,y=>$y,width=>$unit_w,height=>$unit_h,style=>"fill:$_;stroke:#000;stroke-width:1");
        $x += $unit_w;
    }

    $y += $unit_h + $spacing;
}

open OUT,">","ColPaletters.svg";
print OUT $svg->xmlify;
close OUT;

#-------------------------------------------------------------------------------
#  sub func
#-------------------------------------------------------------------------------
sub read_conf {
    my $file = shift;

    open my $fh_conf , $file or die $!;
    while(<$fh_conf>){
        chomp;
        next if (/^#/ || $_ eq "");
        s/#(.+)$//;
        s/\s//g;

        my ($name , $val) = split /=/ , $_;
        $conf{$name} = hex_rgb($val);
    }
    close $fh_conf;

    return %conf;
}

sub fetch_pals {
    my %list = map { s/\-\w+$//; $_ => 1 } @_;

    my @pals = map { $_->[0]  } 
               sort { 
#                      $a->[1] cmp $b->[1] 
#                               || 
#                      $a->[2] <=> $b->[2] } 
                      $a->[2] <=> $b->[2] 
                               || 
                      $a->[1] cmp $b->[1] } 
               map { [ $_ , split /-/ , $_ ] }  keys %list;
    
    return @pals;
}
