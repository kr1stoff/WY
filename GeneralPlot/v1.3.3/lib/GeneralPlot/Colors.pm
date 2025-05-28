package GeneralPlot::Colors;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.04.09 16:28:46          |
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
our @EXPORT = qw(fetch_color);
our @EXPORT_OK = qw();

use lib "$Bin";
use lib "$Bin/..";

use GeneralPlot::Envs;
use GeneralPlot::Debug;

#-------------------------------------------------------------------------------
#  fetch preset color 
#-------------------------------------------------------------------------------

sub fetch_color {

    # a R version written by xushuyang
    # from now on, color sets would update in ../../src/Rlib/Colors/Colors.R

    my %confs = @_;
    !defined $confs{'-type'} and ERROR("fetch_color", "<-type> is not defined. ");
    !defined $confs{'-num'}  and $confs{'-num'} = 0;
    !defined $confs{'-tag'}  and $confs{'-tag'} = "default";
    print join("\t",'fetch_color', $confs{'-type'}, $confs{'-num'}, $confs{'-tag'}) . "\n";

    my $colors = `$SOFT{Rscript} $SOFT{colors_run} $confs{'-type'} $confs{'-tag'} $confs{'-num'}`;
    my @sets = split /,/, $colors;
    return [@sets];

}

sub fetch_color_v3 {
    
    # a update version for more compatibility
    my %confs = @_;
    !defined $confs{'-type'} and ERROR("fetch_color", "<-type> is not defined. ");
    !defined $confs{'-num'}  and $confs{'-num'} = 0;
    !defined $confs{'-tag'}  and $confs{'-tag'} = "default";
    print join("\t",'fetch_color', $confs{'-type'}, $confs{'-num'}, $confs{'-tag'}) . "\n";
    
    my $sets = "color_$confs{'-type'}_$confs{'-tag'}";
    !%{$sets} and ERROR("fetch_color", "<-type $confs{'-type'}> is not right. ");

    if(defined ${$sets}{$confs{'-num'}}) {
        return ${$sets}{$confs{'-num'}};
    }
    else {
        for my $key (sort {scalar(@{$sets->{$a}}) <=> scalar(@{$sets->{$b}})} keys %{$sets}) {
            if(scalar(@{$sets->{$key}}) >= $confs{'-num'}) {
                my @sets = @{$sets->{$key}}[0 .. ($confs{'-num'}-1)];
                return [@sets];
            }
        }

        my $color_str = join(",", @{$sets->{0}});
        my $colors_par = `$SOFT{Rscript} $SOFT{get_color_pal} '$color_str' $confs{'-num'}`;
        my @sets = split /,/, $colors_par;
        return [@sets];
    }

}

sub _fetch_color_v2 {
    
    # a update version for more compatibility
    my %confs = @_;
    !defined $confs{'-type'} and ERROR("fetch_color", "<-type> is not defined. ");
    !defined $confs{'-num'}  and $confs{'-num'} = 0;
    !defined $confs{'-tag'}  and $confs{'-tag'} = "default";
    print join("\t",'fetch_color', $confs{'-type'}, $confs{'-num'}, $confs{'-tag'}) . "\n";
    
    my $sets = "color_$confs{'-type'}_$confs{'-tag'}";
    !%{$sets} and ERROR("fetch_color", "<-type $confs{'-type'}> is not right. ");

    if(defined ${$sets}{$confs{'-num'}}) {
        return ${$sets}{$confs{'-num'}};
    }
    else {
        for my $key (sort {scalar(@{$sets->{$a}}) <=> scalar(@{$sets->{$b}})} keys %{$sets}) {
            if(scalar(@{$sets->{$key}}) >= $confs{'-num'}) {
                return $sets->{$key};
            }
        }
        return $sets->{0};
    }

}

sub _fetch_color_v1 {
    
    my %confs = @_;
    !defined $confs{'-type'} and ERROR("fetch_color", "<-type> is not defined. ");
    !defined $confs{'-num'}  and $confs{'-num'} = 0;
    !defined $confs{'-tag'}  and $confs{'-tag'} = "default";
    
    my $sets = "color_$confs{'-type'}_$confs{'-tag'}";
    !%{$sets} and ERROR("fetch_color", "<-type $confs{'-type'}> is not right. ");

    ${$sets}{$confs{'-num'}} ? 
    ${$sets}{$confs{'-num'}} : 
    ERROR("fetch_color", "<-num $confs{'-num'}> is not defined in default color sets. ");

}

#-------------------------------------------------------------------------------
#  preset color pattern for venn chart 
#-------------------------------------------------------------------------------

our %color_venn_default = (

    2 => [qw/#FFBDC0 #C7D4EE/],
    3 => [qw/#69A4F9 #FFCC66 #FCB4B4/],
    4 => [qw/#FFCC66 #BFE046 #FCB4B4 #69A4F9/],
    5 => [qw/#BFE046 #69A4F9 #ACB9EA #FCB4B4 #FFCC66/],
    6 => [qw/#FCB4B4 #FFCC66 #BFE046 #28C580 #69A4F9 #ACB9EA/],
    7 => [qw/#FCB4B4 #FFCC66 #BFE046 #28C580 #69A4F9 #ACB9EA #C3C3C3/],
    0 => [qw/#FCB4B4 #FFCC66 #BFE046 #28C580 #69A4F9 #ACB9EA #C3C3C3/],

);

#-------------------------------------------------------------------------------
#  preset color pattern for line chart 
#-------------------------------------------------------------------------------
our %color_line_default = (

    1  => [qw/#00468B/],
    2  => [qw/#D32421 #00468B/],
    3  => [qw/#0A7C2E #00468B #D32421/],
    4  => [qw/#0A7C2E #00468B #925E9F #D32421/],
    5  => [qw/#00468B #0A7C2E #0099B4 #925E9F #AD002A/],
    8  => [qw/#00468B #925E9F #0099B4 #0A7C2E #FDAF91 #FD0000 #AD002A #ADB6B6/],
    10 => [qw/#00468B #925E9F #0099B4 #3DB88C #0A7C2E #FDAF91 #FD0000 #AD002A 
              #ADB6B6 #1B1919/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #E67E74 #FF7777 #FD0000 #AD002A #792244 
              #AD556B #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

);

our %color_line_defaultv1 = (

    1 => [qw/#3852A4/],
    2 => [qw/#3852A4 #9C3487/],
    3 => [qw/#3852A4 #9C3487 #55A0FB/],
    4 => [qw/#3852A4 #9C3487 #55A0FB #28C580/],
    0 => [qw/#A00100 #CB2507 #F8403F #FF8180 #FFB3A1 #D7B0B0 #8DAAD3 #55A0FB
             #5758FB #0037FE #00177F #7E369A #9C3487 #477189 #26AAA7 #23AD79
             #028D44 #26711E #979E0C #CDAC10/],

);

#-------------------------------------------------------------------------------
#  preset color pattern for pie chart 
#-------------------------------------------------------------------------------
our %color_pie_default = (

    1  => [qw/#0099B4/],
    2  => [qw/#42B540 #0099B4/],
    3  => [qw/#00468B #0099B4 #42B540/],
    4  => [qw/#00468B #42B540 #0099B4 #EDE447/],
    5  => [qw/#00468B #0099B4 #42B540 #EDE447 #FF7777/],
    8  => [qw/#00468B #0099B4 #76D1B1 #42B540 #EDE447 #FF7777 #AD002A #759EDD/],
    10 => [qw/#ADB6B6 #1B1919 #00468B #0099B4 #76D1B1 #42B540 #EDE447 #FF7777 
              #AD002A #759EDD/],
    20 => [qw/#ADB6B6 #1B1919 #00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 
              #42C1BB #76D1B1 #42B540 #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 
              #FD0000 #AD002A #AE8691 #AE8691/],
    30 => [qw/#ADB6B6 #1B1919 #00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F 
              #759EDD #76C8DC #0099B4 #42C1BB #76D1B1 #0F8074 #28AA6C #42B540 
              #B8D24D #EDE447 #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A 
              #792244 #AD556B #AE8691 #CE9573 #B09F91 #756455/],
    0  => [qw/#c9caca #3b3230 #2e0a4a #7a1b6c #15ad68 #ded531 #db9421 #14b5b5
              #ede893 #76cfed #4599de #db2830 #5F64AC #B271AD #eda4d3 #d33c67
              #EC6925 #a155f9 #70F2D3 #6FBA33 #EDAC2E #096d42 #4ec4b5 #a36924 
              #125fb2 #7350EB #891a3a #bf109a #E8851F #e77def #4ebee5 #69a4f9 
              #f9cfa5 #13D9B1 #bfe046/],
);

our %color_pie_defaultv1 = (

    1 => [qw/#FFCAA8/],
    2 => [qw/#FFCAA8 #80B0DB/],
    3 => [qw/#F27762 #4BDAB5 #4394CD/],
    4 => [qw/#FFB775 #23D3AC #97D9FF #A9B8F1/],
    5 => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8/],
    6 => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8 #FDB462/],
    7 => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8 #FDB462 #B3DE69/],
    8 => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8 #FDB462 #B3DE69 #FCCDE5/],
    0 => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8 #FDB462 #B3DE69 #FCCDE5
             #D9D9D9 #BE82BF #CCEBC5 #FFF673/],
);

our %color_twopie_default = (

    1  => [qw/#2771A7/],
    2  => [qw/#D32421 #2771A7/],
    3  => [qw/#3A9736 #2771A7 #D32421/],
    4  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421/],
    5  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421 #F3BB6F/],
    8  => [qw/#D32421 #F09594 #2771A7 #3A9736 #F3BB6F #C6AFD1 #831D20 #A2C8DC/],
    10 => [qw/#A2C8DC #F09594 #2771A7 #C6AFD1 #D32421 #831D20 #3A9736 #F3BB6F 
              #A3A49E #040000/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B 
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

);

#-------------------------------------------------------------------------------
#  preset color pattern for bar and hist chart 
#-------------------------------------------------------------------------------
our %color_bar_default = (

    1  => [qw/#2771A7/],
    2  => [qw/#D32421 #2771A7/],
    3  => [qw/#3A9736 #2771A7 #D32421/],
    4  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421/],
    5  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421 #F3BB6F/],

    8  => [qw/#00468B #925E9F #0099B4 #0A7C2E #FDAF91 #FD0000 #AD002A #ADB6B6/],
    10 => [qw/#00468B #925E9F #0099B4 #3DB88C #0A7C2E #FDAF91 #FD0000 #AD002A 
              #ADB6B6 #1B1919/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

);

our %color_bar_defaultv1 = (

    1  => [qw/#3852A4/],
    2  => [qw/#3852A4 #9C3487/],
    3  => [qw/#3852A4 #9C3487 #55A0FB/],
    4  => [qw/#3852A4 #9C3487 #55A0FB #28C580/],
    5  => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8/],
    6  => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8 #FDB462/],
    7  => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8 #FDB462 #B3DE69/],
    8  => [qw/#80B1D2 #FB8071 #BEBBDA #FEFFB3 #8DD3C8 #FDB462 #B3DE69 #FCCDE5/],
    0  => [qw/#A00100 #CB2507 #F8403F #FF8180 #FFB3A1 #D7B0B0 #8DAAD3 #55A0FB
              #5758FB #0037FE #00177F #7E369A #9C3487 #477189 #26AAA7 #23AD79
              #028D44 #26711E #979E0C #CDAC10/],

);

our %color_bar_stack = (

    1  => [qw/#42B540/],
    2  => [qw/#42B540 #EDE447/],
    3  => [qw/#42B540 #EDE447 #FF7777/],
    4  => [qw/#00468B #42B540 #EDE447 #FF7777/],
    5  => [qw/#00468B #42B540 #EDE447 #759EDD #FF7777/],
    8  => [qw/#00468B #0099B4 #76D1B1 #42B540 #EDE447 #FF7777 #AD002A #759EDD/],
    10 => [qw/#ADB6B6 #1B1919 #00468B #0099B4 #76D1B1 #42B540 #EDE447 #FF7777 
              #AD002A #759EDD/],
    20 => [qw/#ADB6B6 #1B1919 #00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 
              #42C1BB #76D1B1 #42B540 #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 
              #FD0000 #AD002A #AE8691 #CE9573/],
    30 => [qw/#ADB6B6 #1B1919 #00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F 
              #759EDD #76C8DC #0099B4 #42C1BB #76D1B1 #0F8074 #28AA6C #42B540 
              #B8D24D #EDE447 #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A 
              #792244 #AD556B #AE8691 #CE9573 #B09F91 #756455/],
    0  => [qw/#c9caca #3b3230 #2e0a4a #7a1b6c #15ad68 #ded531 #db9421 #14b5b5
              #ede893 #76cfed #4599de #db2830 #5F64AC #B271AD #eda4d3 #d33c67
              #EC6925 #a155f9 #70F2D3 #6FBA33 #EDAC2E #096d42 #4ec4b5 #a36924 
              #125fb2 #7350EB #891a3a #bf109a #E8851F #e77def #4ebee5 #69a4f9 
              #f9cfa5 #13D9B1 #bfe046/],
);

our %color_bar_stackv1 = (

    1  => [qw/#0099B4/],
    2  => [qw/#42B540 #0099B4/],
    3  => [qw/#00468B #0099B4 #42B540/],
    4  => [qw/#00468B #42B540 #0099B4 #EDE447/],
    5  => [qw/#00468B #0099B4 #42B540 #EDE447 #FF7777/],
    8  => [qw/#00468B #0099B4 #76D1B1 #42B540 #EDE447 #FF7777 #AD002A #759EDD/],
    10 => [qw/#ADB6B6 #1B1919 #00468B #0099B4 #76D1B1 #42B540 #EDE447 #FF7777 
              #AD002A #759EDD/],
    20 => [qw/#ADB6B6 #1B1919 #00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 
              #42C1BB #76D1B1 #42B540 #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 
              #FD0000 #AD002A #AE8691 #AE8691/],
    30 => [qw/#ADB6B6 #1B1919 #00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F 
              #759EDD #76C8DC #0099B4 #42C1BB #76D1B1 #0F8074 #28AA6C #42B540 
              #B8D24D #EDE447 #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A 
              #792244 #AD556B #AE8691 #CE9573 #B09F91 #756455/],
    0  => [qw/#c9caca #3b3230 #2e0a4a #7a1b6c #15ad68 #ded531 #db9421 #14b5b5
              #ede893 #76cfed #4599de #db2830 #5F64AC #B271AD #eda4d3 #d33c67
              #EC6925 #a155f9 #70F2D3 #6FBA33 #EDAC2E #096d42 #4ec4b5 #a36924 
              #125fb2 #7350EB #891a3a #bf109a #E8851F #e77def #4ebee5 #69a4f9 
              #f9cfa5 #13D9B1 #bfe046/],
);

our %color_hist_default = (

    1  => [qw/#00468B/],
    2  => [qw/#D32421 #00468B/],
    3  => [qw/#0A7C2E #00468B #D32421/],
    4  => [qw/#0A7C2E #00468B #925E9F #D32421/],
    5  => [qw/#00468B #0A7C2E #0099B4 #925E9F #AD002A/],
    8  => [qw/#00468B #925E9F #0099B4 #0A7C2E #FDAF91 #FD0000 #AD002A #ADB6B6/],
    10 => [qw/#00468B #925E9F #0099B4 #3DB88C #0A7C2E #FDAF91 #FD0000 #AD002A 
              #ADB6B6 #1B1919/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

);

#-------------------------------------------------------------------------------
#  preset color pattern for heatmap chart 
#-------------------------------------------------------------------------------
our %color_heatmap_default = (

    0 => [qw/#2F70AD #FFFFFF #BA2831/],
    1 => [qw/#4785B6 #FFFFFF #FF1717/],

    2 => [qw/#5B8089 #FFFFFF #009933/],
    3 => [qw/#5B8089 #FDECBE #009933/],
    4 => [qw/#3F99CB #FFFFFD #A1CD46/],
    5 => [qw/#FFFD36 #FDFFFE #2C9B68/],
    6 => [qw/#013F84 #FDECBE #02908B/],

);

our %color_heatmap_corr = (

    0 => [qw/#FFFFFF #C4DEEC #2166AC/],
    1 => [qw/#4785B6 #FFFFFF #FF1717/],

    2 => [qw/#5B8089 #FFFFFF #009933/],
    3 => [qw/#5B8089 #FDECBE #009933/],
    4 => [qw/#3F99CB #FFFFFD #A1CD46/],
    5 => [qw/#FFFD36 #FDFFFE #2C9B68/],
    6 => [qw/#013F84 #FDECBE #02908B/],
    7 => [qw/#2F70AD #FFFFFF #BA2831/],

);

#-------------------------------------------------------------------------------
#  preset color pattern for dot chart 
#-------------------------------------------------------------------------------
our %color_dot_default = (

    1  => [qw/#2771A7/],   
    2  => [qw/#D32421 #2771A7/],
    3  => [qw/#3A9736 #2771A7 #D32421/],
    4  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421/],
    5  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421 #F3BB6F/],
    8  => [qw/#D32421 #F09594 #2771A7 #3A9736 #F3BB6F #C6AFD1 #831D20 #A2C8DC/],
    10 => [qw/#A2C8DC #F09594 #2771A7 #C6AFD1 #D32421 #831D20 #3A9736 #F3BB6F 
              #A3A49E #040000/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B 
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

    volcano => [qw/#F9766D #00000032 #609DFF/],

);

our %color_dot_defaultv1 = (

    1  => [qw/#3852A4/],
    2  => [qw/#3852A4 #9C3487/],
    3  => [qw/#3852A4 #9C3487 #55A0FB/],
    4  => [qw/#3852A4 #9C3487 #55A0FB #28C580/],
    5  => [qw/#3852A4 #9C3487 #55A0FB #28C580 #FFCC00/],
    10 => [qw/#3852A4 #9C3487 #55A0FB #28C580 #FFCC00 #FF0000 #D5F30B #26D7AE
              #97D9FF #E01E84/],
    18 => [qw/#FF0000 #FF7300 #FFAF00 #FFEC00 #D5F30B #52D726 #1BAA2F #2DCB75
              #26D7AE #7CDDDD #5FB7D4 #97D9FF #007ED6 #8399EB #8E6CEF #9C46D0
              #C758D0 #E01E84/],
    0  => [qw/#A00100 #CB2507 #F8403F #FF8180 #FFB3A1 #D7B0B0 #8DAAD3 #55A0FB
              #5758FB #0037FE #00177F #7E369A #9C3487 #477189 #26AAA7 #23AD79
              #028D44 #26711E #979E0C #CDAC10/],
    volcano => [qw/#F9766D #00000032 #609DFF/],

);

#-------------------------------------------------------------------------------
#  preset color pattern for box chart 
#-------------------------------------------------------------------------------
our %color_box_default = (

    1  => [qw/#2771A7/],
    2  => [qw/#D32421 #2771A7/],
    3  => [qw/#3A9736 #2771A7 #D32421/],
    4  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421/],
    5  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421 #F3BB6F/],
    8  => [qw/#D32421 #F09594 #2771A7 #3A9736 #F3BB6F #C6AFD1 #831D20 #A2C8DC/],
    10 => [qw/#A2C8DC #F09594 #2771A7 #C6AFD1 #D32421 #831D20 #3A9736 #F3BB6F 
              #A3A49E #040000/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B 
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

    split  => [qw/#56B4E9 #E69F00/],
    split2 => [qw/#F9766D #609DFF/],

);

our %color_box_defaultv1 = (

    1  => [qw/#3852A4/],
    2  => [qw/#3852A4 #9C3487/],
    3  => [qw/#3852A4 #9C3487 #55A0FB/],
    4  => [qw/#3852A4 #9C3487 #55A0FB #28C580/],
    18 => [qw/#FF0000 #FF7300 #FFAF00 #FFEC00 #D5F30B #52D726 #1BAA2F #2DCB75
              #26D7AE #7CDDDD #5FB7D4 #97D9FF #007ED6 #8399EB #8E6CEF #9C46D0
              #C758D0 #E01E84/],
    20 => [qw/#A00100 #CB2507 #F8403F #FF8180 #FFB3A1 #D7B0B0 #8DAAD3 #55A0FB
              #5758FB #0037FE #00177F #7E369A #9C3487 #477189 #26AAA7 #23AD79
              #028D44 #26711E #979E0C #CDAC10/],
    0  => [qw/#A00100 #CB2507 #F8403F #FF8180 #FFB3A1 #D7B0B0 #8DAAD3 #55A0FB
              #5758FB #0037FE #00177F #7E369A #9C3487 #477189 #26AAA7 #23AD79
              #028D44 #26711E #979E0C #CDAC10/],

    split  => [qw/#56B4E9 #E69F00/],
    split2 => [qw/#F9766D #609DFF/],

);

#-------------------------------------------------------------------------------
#  preset color pattern for all chart 
#-------------------------------------------------------------------------------
our %color_all_default = (

    1  => [qw/#2771A7/],   
    2  => [qw/#D32421 #2771A7/],
    3  => [qw/#3A9736 #2771A7 #D32421/],
    4  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421/],
    5  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421 #F3BB6F/],
    8  => [qw/#D32421 #F09594 #2771A7 #3A9736 #F3BB6F #C6AFD1 #831D20 #A2C8DC/],
    10 => [qw/#A2C8DC #F09594 #2771A7 #C6AFD1 #D32421 #831D20 #3A9736 #F3BB6F 
              #A3A49E #040000/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B 
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

);

our %color_all_defaultv1 = (

    1  => [qw/#0099B4/],
    2  => [qw/#42B540 #0099B4/],
    3  => [qw/#42B540 #0099B4 #925E9F/],
    4  => [qw/#00468B #42B540 #0099B4 #925E9F/],
    5  => [qw/#00468B #42B540 #0099B4 #925E9F #AD002A/],
    8  => [qw/#00468B #925E9F #0099B4 #42B540 #FDAF91 #FD0000 #AD002A #ADB6B6/],
    10 => [qw/#00468B #925E9F #0099B4 #3DB88C #42B540 #FDAF91 #FD0000 #AD002A
              #ADB6B6 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #0099B4 #925E9F #759EDD #0099B4 #42C1BB
              #76D1B1 #42B540 #B8D24D #EDE447 #FAB158 #FF7777 #FD0000 #AD002A
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC
              #0099B4 #42C1BB #76D1B1 #0F8074 #28AA6C #42B540 #B8D24D #EDE447
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #0099B4 #925E9F #759EDD #0099B4 #42C1BB
              #76D1B1 #42B540 #B8D24D #EDE447 #FAB158 #FF7777 #FD0000 #AD002A
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    
);

# dark color set
our %color_all_dark = (

    1  => [qw/#0A7C2E/],
    2  => [qw/#D32421 #00468B/],
    3  => [qw/#0A7C2E #00468B #D32421/],
    4  => [qw/#0A7C2E #00468B #925E9F #D32421/],
    5  => [qw/#00468B #0A7C2E #0099B4 #925E9F #AD002A/],
    8  => [qw/#00468B #925E9F #0099B4 #0A7C2E #FDAF91 #FD0000 #AD002A #ADB6B6/],
    10 => [qw/#00468B #925E9F #0099B4 #3DB88C #0A7C2E #FDAF91 #FD0000 #AD002A 
              #ADB6B6 #1B1919/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B 
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

);

# light color set
our %color_all_light = (

    1  => [qw/#2771A7/],    
    2  => [qw/#D32421 #2771A7/],
    3  => [qw/#3A9736 #2771A7 #D32421/],
    4  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421/],
    5  => [qw/#3A9736 #2771A7 #C6AFD1 #D32421 #F3BB6F/],
    8  => [qw/#D32421 #F09594 #2771A7 #3A9736 #F3BB6F #C6AFD1 #831D20 #A2C8DC/],
    10 => [qw/#A2C8DC #F09594 #2771A7 #C6AFD1 #D32421 #831D20 #3A9736 #F3BB6F 
              #A3A49E #040000/],
    15 => [qw/#00468B #925E9F #759EDD #0099B4 #0A7C2E #B8D24D #EDE447 #FAB158 
              #FF7777 #FD0000 #AD002A #AE8691 #ADB6B6 #CE9573 #1B1919/],
    20 => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],
    30 => [qw/#00468B #5377A7 #3B81AB #5C298F #6C6DA4 #925E9F #759EDD #76C8DC 
              #0099B4 #42C1BB #76D1B1 #0F8074 #0A7C2E #28AA6C #B8D24D #EDE447 
              #FAB158 #FDAF91 #E67E74 #FF7777 #FD0000 #AD002A #792244 #AD556B 
              #AE8691 #CE9573 #B09F91 #ADB6B6 #4C4E4E #1B1919/],
    0  => [qw/#00468B #5377A7 #6C6DA4 #925E9F #759EDD #0099B4 #42C1BB #76D1B1 
              #0A7C2E #B8D24D #EDE447 #FAB158 #FDAF91 #FF7777 #FD0000 #AD002A 
              #AE8691 #ADB6B6 #4C4E4E #1B1919/],

);

#-------------------------------------------------------------------------------
#  preset color pattern from ggsci 
#-------------------------------------------------------------------------------
our %color_ggsci_npg = (

    0  => [qw/#E64B35 #4DBBD5 #00A087 #3C5488 #F39B7F #8491B4 #91D1C2 #DC0000 
              #7E6148 #B09C85/],
);

our %color_ggsci_lancet = (

    0  => [qw/#00468B #ED0000 #42B540 #0099B4 #925E9F #FDAF91 #AD002A #ADB6B6 
              #1B1919/],

);

return 1;

