package SBV::STONE::SYMBOL;
#------------------------------------------------+
#    [APM] This moudle is generated by amp.pl    |
#    [APM] Creat time: 2013-06-03 15:49:46       |
#------------------------------------------------+
=pod

=head1 Name

SBV::STONE::SYMBOL

=head1 Synopsis

This module is not meant to be used directly

=head1 Feedback

Author: Peng Ai
Email:  aipeng0520@163.com

=head1 Version

Version history

=head2 v1.0

Date: 2013-06-03 15:49:46

=cut

use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
#our @EXPORT    = qw(symbol);

use Math::Round;

use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/../";
use lib "$FindBin::RealBin/../lib";

use SBV;
use SBV::DEBUG;
use SBV::Constants;

=pod

=head1 Name

new -- create the new symbol defs in svg

=head1 Run

new($pch,%par)

=head1 Parameter

$pch -- the shape of the symbol

width -- the symbol width 
height -- the symbol height
fill -- the fill color 
color -- the lines/stroke color

=head1 PCH

0 -- rect
1 -- circle

=cut
sub new
{
    my ($pch,%par) = @_;
    
    # init the SYMBOL paramter
    my ($w,$h,$size,$fill,$color,$swidth,%newpar) = doInitPar(%par);
    my $r = $w > $h ? $h/2 : $w/2;
    
    my $id = "symbol$SBV::idnum";
    $SBV::idnum ++;
    my $parent = $par{parent} || $SBV::defs;
    my $symbol = $parent->symbol(id=>$id,class=>"symbol",viewBox=>"0 0 $w $h");
    
    # styles
    my $style1 = {stroke=>$color,fill=>$fill,'stroke-width'=>$swidth};
    my $style2 = {stroke=>$color,'stroke-width'=>$swidth};
    my $style3 = {stroke=>$color,fill=>"none",'stroke-width'=>$swidth};
    my $style4 = {fill=>$color,'stroke-width'=>0}; # for points on the lines 
    $style1->{'fill-opacity'} = $par{opacity} if (exists $par{opacity});

    # specific points coordinate
    my $cx = $w/2;    my $cy = $h/2; # central point
    
    my $ox = $swidth;    my $oy = $swidth; # top left point
    my $x2 = $w - $swidth; my $y2 = $h - $swidth; # bottom right point
    
    #my $ox = 0;    my $oy = 0; # top left point
    #my $x2 = $w; my $y2 = $h; # bottom right point
    
    $w -= 2*$ox;
    $h -= 2*$oy;
    $r -= $ox > $oy ? $oy : $ox;
    
    if (35 != $pch){
        $size = 1 if ($size > 1);
        $ox += (1-$size)*$w;
        $oy += (1-$size)*$h;
        $x2 -= (1-$size)*$w;
        $y2 -= (1-$size)*$h;
        $r *= $size;
        $w *= $size;
        $h *= $size;
    }else{
        $size = 1/3 if ($size == 1)
    }
    
    # draw background
    #my $symbol_conf = {ox=>0,oy=>$h,oty=>0,tw=>$w,th=>$h,
    #    background=>$newpar{background},border=>$newpar{border}};
    #SBV::DRAW::background($symbol_conf,$symbol);
    
    # draw symbol
    if (0 == $pch) # rect
    {
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,style=>$style1);
    }
    elsif (1 == $pch) # circle
    {
        $symbol->circle(cx=>$cx,cy=>$cy,r=>$r,style=>$style1);
    }
    elsif (2 == $pch) # diomand
    {
        $symbol->path(d=>"M$cx $oy L$x2 $cy L$cx $y2 L$ox $cy Z",style=>$style1);
    }
    elsif (3 == $pch) # regular triangle (the oriental is top)
    {
        my ($x1,$y1,$x2,$y2,$x3,$y3);
        $x1 = $cx; $y1 = $cy - $r;
        $x2 = nearest 0.0001 , $cx - $r*$SQRT3/2; $y2 = $cy + $r/2;
        $x3 = nearest 0.0001 , $cx + $r*$SQRT3/2; $y3 = $cy + $r/2;
        $symbol->path(d=>"M$x1 $y1 L$x2 $y2 L$x3 $y3 Z",style=>$style1);
    }
    elsif (4 == $pch) # regular triangle (the oriental is bottom)
    {
        my ($x1,$y1,$x2,$y2,$x3,$y3);
        $x1 = $cx; $y1 = $cy + $r;
        $x2 = nearest 0.0001 , $cx - $r*$SQRT3/2; $y2 = $cy - $r/2;
        $x3 = nearest 0.0001 , $cx + $r*$SQRT3/2; $y3 = $cy - $r/2;
        $symbol->path(d=>"M$x1 $y1 L$x2 $y2 L$x3 $y3 Z",style=>$style1);
    }
    elsif (5 == $pch) # cross mark: '+'
    {
        $symbol->line(x1=>$cx,y1=>$ox,x2=>$cx,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,y1=>$cy,x2=>$x2,y2=>$cy,style=>$style2);
    }
    elsif (6 == $pch) # error mark: 'X'
    {
        $symbol->line(x1=>$ox,y1=>$oy,x2=>$x2,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,y1=>$y2,x2=>$x2,y2=>$oy,style=>$style2);
    }
    elsif (7 == $pch) # '+' and 'X'
    {
        # '+'
        $symbol->line(x1=>$cx,y1=>$ox,x2=>$cx,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,y1=>$cy,x2=>$x2,y2=>$cy,style=>$style2);
        
        # 'X'
        $symbol->line(x1=>$ox,y1=>$oy,x2=>$x2,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,y1=>$y2,x2=>$x2,y2=>$oy,style=>$style2);
    }
    elsif (8 == $pch) # rect + '+'
    {
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,style=>$style3);
        # '+'
        $symbol->line(x1=>$cx,y1=>$ox,x2=>$cx,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,y1=>$cy,x2=>$x2,y2=>$cy,style=>$style2);
    }
    elsif (9 == $pch) # rect and 'X'
    {
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,style=>$style3);
        
        # 'X'
        $symbol->line(x1=>$ox,y1=>$oy,x2=>$x2,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,y1=>$y2,x2=>$x2,y2=>$oy,style=>$style2);
    }
    elsif (10 == $pch) # rect + 'V'
    {
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,style=>$style3);
        # 'V'
        $symbol->line(x1=>$ox,y1=>$oy,x2=>$cx,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$cx,y1=>$y2,x2=>$x2,y2=>$oy,style=>$style2);
    }
    elsif (11 == $pch) # circle and '+'
    {
        $symbol->circle(cx=>$cx,cy=>$cy,r=>$r,style=>$style3);    
        # '+'
        $symbol->line(x1=>$cx,y1=>$cy-$r,x2=>$cx,y2=>$cy+$r,style=>$style2);
        $symbol->line(x1=>$cx-$r,y1=>$cy,x2=>$cx+$r,y2=>$cy,style=>$style2);
    }
    elsif (12 == $pch) # circle + 'X'
    {
        $symbol->circle(cx=>$cx,cy=>$cy,r=>$r,style=>$style3);    
        # 'X'
        $symbol->line(x1=>$cx-$r,y1=>$cy-$r,x2=>$cx+$r,y2=>$cy+$r,style=>$style2);
        $symbol->line(x1=>$cx-$r,y1=>$cy+$r,x2=>$cx+$r,y2=>$cy-$r,style=>$style2);
    }
    elsif (13 == $pch) # diamond and '+'
    {
        $symbol->path(d=>"M$cx $oy L$x2 $cy L$cx $y2 L$ox $cy Z",style=>$style3);
        
        # '+'
        $symbol->line(x1=>$cx,y1=>$ox,x2=>$cx,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,y1=>$cy,x2=>$x2,y2=>$cy,style=>$style2);
    }
    elsif (14 == $pch) # 2 + 6 (two regular triangle)
    {
        # top triangle
        my ($x1,$y1,$x2,$y2,$x3,$y3);
        $x1 = $cx; $y1 = $cy - $r;
        $x2 = nearest 0.0001 , $cx - $r*$SQRT3/2; $y2 = $cy + $r/2;
        $x3 = nearest 0.0001 , $cx + $r*$SQRT3/2; $y3 = $cy + $r/2;
        $symbol->path(d=>"M$x1 $y1 L$x2 $y2 L$x3 $y3 Z",style=>$style1);
        
        # bottom triangle
        my ($sx1,$sy1,$sx2,$sy2,$sx3,$sy3);
        $sx1 = $cx; $sy1 = $cy + $r;
        $sx2 = nearest 0.0001 , $cx - $r*$SQRT3/2; $sy2 = $cy - $r/2;
        $sx3 = nearest 0.0001 , $cx + $r*$SQRT3/2; $sy3 = $cy - $r/2;
        $symbol->path(d=>"M$sx1 $sy1 L$sx2 $sy2 L$sx3 $sy3 Z",style=>$style1);
        $symbol->path(d=>"M$x1 $y1 L$x2 $y2 L$x3 $y3 Z",style=>$style3);
        $symbol->path(d=>"M$sx1 $sy1 L$sx2 $sy2 L$sx3 $sy3 Z",style=>$style3);

    }
    elsif (15 == $pch) # five-pointed star 
    {
        fpstar($r,$cx,$cy,$style1,$symbol);
    }
    elsif (16 == $pch) #six-pointed star
    {
        spstar($r,$cx,$cy,$style1,$symbol);    
    }
    elsif (17 == $pch) #line
    {
        $symbol->line(x1=>$ox,x2=>$x2,y1=>($oy+$y2)/2,y2=>($oy+$y2)/2,style=>$style2);    
    }
    elsif (18 == $pch) # line and a points on it
    {
        $symbol->line(x1=>$ox,x2=>$x2,y1=>($oy+$y2)/2,y2=>($oy+$y2)/2,style=>$style2); 
        $symbol->circle(cx=>$cx,cy=>$cy,r=>3,style=>$style4);
    }
    elsif (19 == $pch) # boxplot
    {
        my $recty = nearest 0.001 , $ox + 0.25*$h;
        my $recty2 = nearest 0.001 , $ox + 0.75*$h;
        $symbol->rect(x=>$ox,y=>$recty,width=>$w,height=>$h/2,style=>$style1);
        $symbol->line(x1=>$ox,x2=>$x2,y1=>$oy+$h/2,y2=>$oy+$h/2,style=>$style2);
        $symbol->line(x1=>$ox+$w/2,x2=>$ox+$w/2,y1=>$oy,y2=>$recty,style=>$style2);
        $symbol->line(x1=>$ox+$w/2,x2=>$ox+$w/2,y1=>$recty2,y2=>$y2,style=>$style2);
    }
    elsif (20 == $pch) # linerange
    {
        $symbol->line(x1=>$cx,x2=>$cx,y1=>$oy,y2=>$y2,style=>$style2);
    }
    elsif (21 == $pch) # pointrange
    {
        $symbol->line(x1=>$cx,x2=>$cx,y1=>$oy,y2=>$y2,style=>$style2);
        $symbol->circle(cx=>$cx,cy=>$cy,r=>3,style=>$style4);
    }
    elsif (22 == $pch) #crossbar
    {
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,style=>$style1);
        $symbol->line(x1=>$ox,x2=>$x2,y1=>$cy,y2=>$cy,style=>$style2);
    }
    elsif (23 == $pch) # errbar
    {
        $symbol->line(x1=>$cx,x2=>$cx,y1=>$oy,y2=>$y2,style=>$style2);
        $symbol->line(x1=>$ox,x2=>$x2,y1=>$oy,y2=>$oy,style=>$style2);
        $symbol->line(x1=>$ox,x2=>$x2,y1=>$y2,y2=>$y2,style=>$style2);
    }
    elsif (24 == $pch) # left pointing pencil
    {
        SBV::DRAW::pencil($ox,$oy,$w,$h,parent=>$symbol,style=>$style1,strand=>-1,arrow_width=>$newpar{arrow_width});
    }
    elsif (25 == $pch) # right pointing pencil
    {
        SBV::DRAW::pencil($ox,$oy,$w,$h,parent=>$symbol,style=>$style1,strand=>1,arrow_width=>$newpar{arrow_width});
    }
    elsif (26 == $pch) # right pointing triangle
    {
        my $px = [$ox,$x2,$ox];
        my $py = [$oy,$cy,$y2];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (27 == $pch) # left pointing triangle
    {
        my $px = [$ox,$x2,$x2];
        my $py = [$cy,$oy,$y2];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (28 == $pch) # right pointing pentagram
    {
        my $px = [$ox,$ox+$w*2/3,$x2,$ox+$w*2/3,$ox];
        my $py = [$oy+$h/4,$oy,$oy+$h/2,$y2,$oy+$h*3/4];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (29 == $pch) # left pointing pentagram
    {
        my $px = [$ox,$ox+$w/3,$x2,$x2,$ox+$w/3];
        my $py = [$oy+$h/2,$oy,$oy+$h/4,$oy+$h*3/4,$y2];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (30 == $pch) # up pointing pentagram
    {
        my $px = [$ox,$cx,$x2,$x2-$w/4,$ox+$w/4];
        my $py = [$oy+$h/3,$oy,$oy+$h/3,$y2,$y2];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (31 == $pch) # down pointing pentagram
    {
        my $px = [$ox+$w/4,$x2-$w/4,$x2,$cx,$ox];
        my $py = [$oy,$oy,$y2-$h/3,$y2,$y2-$h/3];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (32 == $pch) # horizontal hexagon
    {
        my $px = [$ox,$ox+$w/4,$x2-$w/4,$x2,$x2-$w/4,$ox+$w/4];
        my $py = [$cy,$oy,$oy,$cy,$y2,$y2];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (33 == $pch) # vertical hexagon
    {
        my $px = [$ox,$ox,$cx,$x2,$x2,$cx];
        my $py = [$y2-$h/4,$oy+$h/4,$oy,$oy+$h/4,$y2-$h/4,$y2];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (34 == $pch) # octagon
    {
        my $px = [$ox+$w/4,$x2-$w/4,$x2,$x2,$x2-$w/4,$ox+$w/4,$ox,$ox];
        my $py = [$oy,$oy,$oy+$h/4,$y2-$h/4,$y2,$y2,$y2-$h/4,$oy+$h/4];
        polygon($px,$py,$style1,$symbol);
    }
    elsif (35 == $pch) # 1/3 size rect, for GAP in sequence
    {
        my $newh = $size * $h;
        $oy = $oy+$h/2-$newh/2;
        $h = $newh;
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,style=>$style1);
    }
    elsif (36 == $pch) # ellipse
    {
        $symbol->ellipse(cx=>$cx,cy=>$cy,rx=>$w/2,ry=>$h/2,style=>$style1);
    }
    elsif (37 == $pch) # rounded rectangle
    {
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,rx=>$r,ry=>$r,style=>$style1);
    }
    elsif (38 == $pch) # conventional heart
    {
        my $matrix = "";
        my $pathd = "";
        $symbol->path(transform=>$matrix,style=>$style1,d=>$pathd)
    }
    elsif (39 == $pch) # implicit heart
    {
        my $matrix = "";
        my $pathd = "";
        $symbol->path(transform=>$matrix,style=>$style1,d=>$pathd)
    }
    elsif (40 == $pch) # left pointing arrow
    {
        SBV::DRAW::arrow2($ox,$oy,$w,$h,parent=>$symbol,style=>$style1,strand=>-1,
            ratio=>$newpar{ratio},arrow_width=>$newpar{arrow_width});
    }
    elsif (41 == $pch) # right pointing arrow
    {
        SBV::DRAW::arrow2($ox,$oy,$w,$h,parent=>$symbol,style=>$style1,strand=>1,
            ratio=>$newpar{ratio},arrow_width=>$newpar{arrow_width});
    }
    elsif (42 == $pch) # fan (Center on the left)
    {
        my $polar  = SBV::Coordinate::POLAR->new(-$w-2*$ox,$h/2+$oy,parent=>$symbol);
        my $theta  = atan2($h/2+$oy,$w*2+4*$ox);
        my $angle1 = $PI/2 - $theta;
        my $angle2 = $PI/2 + $theta;
        $polar->fan($w+$ox*3,$angle1,$w*2+3*$ox,$angle2,style=>$style1,theta_type=>"theta");
    }
    elsif (43 == $pch) # semicircle 
    {
        my $re_def_r = $w/2 > $h ? $h-$oy*2 : $w/2 - $ox*2;
        my $pathd = "M$ox $y2 A$re_def_r $re_def_r 0 1 1 $x2 $y2";
        $symbol->path(style=>$style1,d=>$pathd)
    }
    elsif (64 == $pch) # 'a' label
    {
        my $font = SBV::Font->new({
            'font-style'  => $par{font_style},
            'font-weight' => $par{font_weight},
            'font-size'   => $par{font_size},
            'font-family' => $par{font_family},
            fill          => $fill
        });
        my $labelH = $font->fetch_text_height;
        my $labelW = $font->fetch_text_width('a');
        my $textx = $ox + $w/2 - $labelW/2;
        my $texty = $oy + $h/2 + $labelH/2;
        $symbol->rect(x=>$ox,y=>$oy,width=>$w,height=>$h,style=>"fill:#ddd;stroke-width:0");
        $symbol->text(x=>$textx,y=>$texty,style=>$font->toStyle)->cdata('a');
    }
    else
    {
        ERROR("no_pch_err",$pch);
    }
    
    if ($newpar{add_line}){
        $symbol->line(x1=>$ox,y1=>$oy+$h,x2=>$ox+$w,y2=>$oy,style=>"stroke:#000;stroke-width:$newpar{add_line}");
    }
    

    if ($newpar{add_text}){
        my $font = SBV::Font->new($newpar{add_text_theme});
        my $texth = $font->fetch_text_height;
        my $text_style = $font->toStyle . "text-anchor:middle;";
        $symbol->text(x=>$ox+$w/2,y=>$oy+$h/2+$texth/2, style=>$text_style)->cdata($newpar{add_text});
    }

    return $id;
}

# five pointed star
sub fpstar
{
    my ($R,$cx,$cy,$style,$father) = @_;

    my @px;
    my @py;

    # 5 external points coordinate
    map {
        my $index = 2*$_;
        $px[$index] = nearest 0.0001 , $cx + sin($index*$TWOPI/10) * $R;
        $py[$index] = nearest 0.0001 , $cy - cos($index*$TWOPI/10) * $R;
    } 0 .. 4;

    # defined the $r for five pointed star internal raduis
    my $r = nearest 0.0001 , $R * sin($TWOPI/20)/cos($TWOPI/10);

    # 5 internal points coordinate
    map {
        my $index = 2*$_ + 1;
        $px[$index] = nearest 0.0001 , $cx + sin($index*$TWOPI/10) * $r;
        $py[$index] = nearest 0.0001 , $cy - cos($index*$TWOPI/10) * $r;
    } 0 .. 4;

    my $points = $father->get_path(x=>\@px,y=>\@py,-type=>'polygon');
    $father->polygon(%$points,style=>$style);
}

# six pointed star
sub spstar
{
    my ($R,$cx,$cy,$style,$father) = @_;
    my @px;
    my @py;

    # 5 external points coordinate
    map { 
        my $index = 2*$_;
        $px[$index] = nearest 0.0001 , $cx + sin($index*$TWOPI/12) * $R;
        $py[$index] = nearest 0.0001 , $cy - cos($index*$TWOPI/12) * $R;
    } 0 .. 5;

    # defined the $r for six pointed star internal raduis
    my $r = $R * sin($TWOPI/12)/sin($TWOPI/6);

    # 5 internal points coordinate
    map {
        my $index = 2*$_ + 1;
        $px[$index] = nearest 0.0001 , $cx + sin($index*$TWOPI/12) * $r;
        $py[$index] = nearest 0.0001 , $cy - cos($index*$TWOPI/12) * $r;
    } 0 .. 5;

    my $points = $father->get_path(x=>\@px,y=>\@py,-type=>'polygon');
    $father->polygon(%$points,style=>$style);
}

sub polygon
{
    my ($px,$py,$style,$parent) = @_;
    my $points = $parent->get_path(x=>$px,y=>$py,-type=>"polygon");
    $parent->polygon(%$points,style=>$style);
}

sub doInitPar
{
    my %par = @_;
    my ($w,$h,$size,$fill,$color,$swidth);
    
    $w = $par{"width"} || 10;
    $h = $par{"height"} || 10;
    $size = $par{"size"} || 1;
    $fill = $par{"fill"} || "none";
    $color = $par{color} ? $par{color} : $par{stroke} ? $par{stroke} : "#000";
    $swidth = exists $par{"stroke_width"} ? $par{"stroke_width"} : 1;
    
    $color = SBV::Colors::fetch_color($color);
    $fill = SBV::Colors::fetch_color($fill);

    $w = nearest 0.01 , $w;
    $h = nearest 0.01 , $h;
    
    $par{background} = "none" unless ($par{background});
    $par{border} = "1111" unless $par{border};
    
    $par{font_style} = "normal" unless $par{font_style};
    $par{font_size} = 16 unless $par{font_size};
    $par{font_weight} = "normal" unless $par{font_weight};
    $par{font_family} = "Arial" unless $par{font_family};

    # just for arrow shape (pch = 40,41)
    $par{ratio} = 0.5 unless $par{ratio};

    # just for arrow shape (pch = 40,41,24,25)
    # the ratio of the arrow height as the arrow width
    $par{arrow_width} = 1 unless $par{arrow_width};

    return ($w,$h,$size,$fill,$color,$swidth,%par);
}

