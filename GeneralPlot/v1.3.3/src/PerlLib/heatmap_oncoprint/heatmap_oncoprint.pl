#!/Bio/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Data::Dumper;
use YAML;
use YAML::Syck;
# Program:       heatmap_oncoprint.pl
# Author:        Liu yubin
# Date:          Wed 28 Oct 2020 10:30:26 AM CST
# Description:   draw a heatmap with ComplexHeatmap oncoPrint
#-----------------------------------------------------------------------------
my $VERSION = "1.1";
my ($conf_file, $matrix_file, $group_file, $group_type, $out_prefix);
my ($rm_empty, $rm_empty_cols, $rm_empty_rows);
my ($show_bar, $show_row_bar, $show_col_bar);
my ($width, $height, $transpose, $title_size, $name_size);
my ($color_list, $color_bg, $color_border);
my ($row_title, $row_title_side, $row_title_size, $row_title_font); 
my ($col_title, $col_title_side, $col_title_size, $col_title_font);
my ($legend_title, $legend_title_size, $legend_title_font, 
    $legend_title_pos, $legend_direction, $legend_side);
my ($show_pct, $pct_size, $pct_font, $pct_digits, $pct_side);
my ($show_row_names, $row_names_size, $row_names_font, 
    $row_names_side, $row_names_rot, $keep_row_order); 
my ($show_col_names, $col_names_size, $col_names_font, 
    $col_names_side, $col_names_rot, $keep_col_order);
my ($Rscript, $convert, $clean, $limit, $help);
my $default_confF = "$Bin/default.yml";

GetOptions(

    "c|conf:s"            => \$conf_file,
    "m|matrix:s"          => \$matrix_file,
    "g|group:s"           => \$group_file,
    "o|outpre:s"          => \$out_prefix,
    
    "l|color_list:s"      => \$color_list,
    "color_bg:s"          => \$color_bg,
    "color_border:s"      => \$color_border,

    "width:s"             => \$width,
    "height:s"            => \$height,
    "transpose:s"         => \$transpose,
    "group_type:s"        => \$group_type,
    "rm_empty:s"          => \$rm_empty,
    "rm_empty_cols:s"     => \$rm_empty_cols,
    "rm_empty_rows:s"     => \$rm_empty_rows,
    "show_bar:s"          => \$show_bar,
    "show_row_bar:s"      => \$show_row_bar,
    "show_col_bar:s"      => \$show_col_bar,
    "title_size:s"        => \$title_size,
    "name_size:s"         => \$name_size,

    "row_title:s"         => \$row_title,
    "row_title_side:s"    => \$row_title_side,
    "row_title_size:s"    => \$row_title_size,
    "row_title_font:s"    => \$row_title_font,
    "col_title:s"         => \$col_title,
    "col_title_side:s"    => \$col_title_side,
    "col_title_size:s"    => \$col_title_size,
    "col_title_font:s"    => \$col_title_font,
    "legend_title:s"      => \$legend_title,
    "legend_title_size:s" => \$legend_title_size,
    "legend_title_font:s" => \$legend_title_font,
    "legend_title_pos:s"  => \$legend_title_pos,
    "legend_direction:s"  => \$legend_direction,
    "legend_side:s"       => \$legend_side,

    "show_pct:s"          => \$show_pct,
    "pct_size:s"          => \$pct_size,
    "pct_font:s"          => \$pct_font,
    "pct_digits:s"        => \$pct_digits,
    "pct_side:s"          => \$pct_side,
    "show_row_names:s"    => \$show_row_names,
    "row_names_size:s"    => \$row_names_size,
    "row_names_font:s"    => \$row_names_font,
    "row_names_side:s"    => \$row_names_side,
    "row_names_rot:s"     => \$row_names_rot,
    "keep_row_order:s"    => \$keep_row_order,
    "show_col_names:s"    => \$show_col_names,
    "col_names_size:s"    => \$col_names_size,
    "col_names_font:s"    => \$col_names_font,
    "col_names_side:s"    => \$col_names_side,
    "col_names_rot:s"     => \$col_names_rot,
    "keep_col_order:s"    => \$keep_col_order,

    "Rscript:s"           => \$Rscript,
    "convert:s"           => \$convert,
    "clean:s"             => \$clean,
    "limit:s"             => \$limit,
    "help!"               => \$help,

);

#-----------------------------------------------------------------------------
if($help) {
    usage();
    exit 0;
}
elsif((!$conf_file) and (!$matrix_file or !$group_file or !$out_prefix)) {
    Warn('no conf file or no input file');
    usage();
    exit 1;
}

#-----------------------------------------------------------------------------
my $conf = {};
my $default_conf = load_conf($default_confF);
if($conf_file and -s $conf_file) {
    $conf = load_conf($conf_file);
}
$conf = init_conf($default_conf, %$conf);

$conf->{matrix} = defined $matrix_file ? $matrix_file : $conf->{matrix};
$conf->{group}  = defined $group_file  ? $group_file  : $conf->{group};
$conf->{outpre} = defined $out_prefix  ? $out_prefix  : $conf->{outpre};

$conf->{width}  = $width ? $width : $conf->{width};
$conf->{height} = $height ? $height : $conf->{height};
$conf->{transpose}  = $transpose ? $transpose : $conf->{transpose};
$conf->{group_type} = $group_type ? $group_type : $conf->{group_type};
$conf->{rm_empty}      = $rm_empty      ? $rm_empty      : $conf->{rm_empty};
$conf->{rm_empty_cols} = $rm_empty_cols ? $rm_empty_cols : $conf->{rm_empty_cols};
$conf->{rm_empty_rows} = $rm_empty_rows ? $rm_empty_rows : $conf->{rm_empty_rows};
$conf->{show_bar}     = $show_bar     ? $show_bar     : $conf->{show_bar};
$conf->{show_row_bar} = $show_row_bar ? $show_row_bar : $conf->{show_row_bar};
$conf->{show_col_bar} = $show_col_bar ? $show_col_bar : $conf->{show_col_bar};
$conf->{title_size} = $title_size ? $title_size : $conf->{title_size};
$conf->{name_size}  = $name_size  ? $name_size  : $conf->{name_size};
$conf->{color_list}   = $color_list ? $color_list : $conf->{color_list};
$conf->{color_bg}     = $color_bg ? $color_bg : $conf->{color_bg};
$conf->{color_border} = $color_border ? $color_border : $conf->{color_border};
$conf->{row_title}      = $row_title ? $row_title : $conf->{row_title};
$conf->{row_title_side} = $row_title_side ? $row_title_side : $conf->{row_title_side};
$conf->{row_title_size} = $row_title_size ? $row_title_size : $conf->{row_title_size};
$conf->{row_title_font} = $row_title_font ? $row_title_font : $conf->{row_title_font};
$conf->{col_title}      = $col_title ? $col_title : $conf->{col_title};
$conf->{col_title_side} = $col_title_side ? $col_title_side : $conf->{col_title_side};
$conf->{col_title_size} = $col_title_size ? $col_title_size : $conf->{col_title_size};
$conf->{col_title_font} = $col_title_font ? $col_title_font : $conf->{col_title_font};
$conf->{legend_title}      = $legend_title ? $legend_title : $conf->{legend_title};
$conf->{legend_title_size} = $legend_title_size ? $legend_title_size : $conf->{legend_title_size};
$conf->{legend_title_font} = $legend_title_font ? $legend_title_font : $conf->{legend_title_font};
$conf->{legend_title_pos}  = $legend_title_pos ? $legend_title_pos : $conf->{legend_title_pos};
$conf->{legend_direction}  = $legend_direction ? $legend_direction : $conf->{legend_direction};
$conf->{legend_side}       = $legend_side ? $legend_side : $conf->{legend_side};
$conf->{show_pct}   = $show_pct ? $show_pct : $conf->{show_pct};
$conf->{pct_size}   = $pct_size ? $pct_size : $conf->{pct_size};
$conf->{pct_font}   = $pct_font ? $pct_font : $conf->{pct_font};
$conf->{pct_digits} = $pct_digits ? $pct_digits : $conf->{pct_digits};
$conf->{pct_side}   = $pct_side ? $pct_side : $conf->{pct_side};
$conf->{show_row_names} = $show_row_names ? $show_row_names : $conf->{show_row_names};
$conf->{row_names_size} = $row_names_size ? $row_names_size : $conf->{row_names_size};
$conf->{row_names_font} = $row_names_font ? $row_names_font : $conf->{row_names_font};
$conf->{row_names_side} = $row_names_side ? $row_names_side : $conf->{row_names_side};
$conf->{row_names_rot}  = $row_names_rot ? $row_names_rot : $conf->{row_names_rot};
$conf->{keep_row_order}  = $keep_row_order ? $keep_row_order : $conf->{keep_row_order};
$conf->{show_col_names} = $show_col_names ? $show_col_names : $conf->{show_col_names};
$conf->{col_names_size} = $col_names_size ? $col_names_size : $conf->{col_names_size};
$conf->{col_names_font} = $col_names_font ? $col_names_font : $conf->{col_names_font};
$conf->{col_names_side} = $col_names_side ? $col_names_side : $conf->{col_names_side};
$conf->{col_names_rot}  = $col_names_rot ? $col_names_rot : $conf->{col_names_rot};
$conf->{keep_col_order}  = $keep_col_order ? $keep_col_order : $conf->{keep_col_order};

$conf->{clean} = $clean ? $clean : $conf->{clean};
$conf->{limit} = $limit ? $limit : $conf->{limit};
$Rscript ||= "/Bio/bin/Rscript-3.6.0";
$convert ||= "/Bio/usr/bin/convert";
#$Rscript ||= "/state/partition1/workspace/software/R-3.6.0/bin/Rscript";
#$convert ||= "/usr/bin/convert";
$conf = check_conf($conf);

#-----------------------------------------------------------------------------
my $outdir  = dirname($conf->{outpre});
my $outname = basename($conf->{outpre});
`mkdir -p $outdir` unless(-d $outdir);

my $matrix = read_matrix($conf->{matrix}, $conf);
my $group  = read_group($conf->{group});
$conf->{group_type} = normalize_matrix($matrix, $group, "$conf->{outpre}.fmt.matrix.txt", $conf->{group_type});
draw_heatmap($matrix, $group, $conf, "$conf->{outpre}.fmt.matrix.txt", $conf->{outpre});
clean_tmp_file($conf);

#-----------------------------------------------------------------------------
sub draw_heatmap {

    my $matrix  = shift;
    my $group   = shift;
    my $conf    = shift;
    my $matrixF = shift;
    my $outpre  = shift;
    my $outdir  = dirname($outpre);
    my $outname = basename($outpre);

    # empty data
    my ($rm_empty_cols, $rm_empty_rows);
    if($conf->{rm_empty_cols} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        $rm_empty_cols = "T";
    } else {
        $rm_empty_cols = "F";
    }
    if($conf->{rm_empty_rows} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        $rm_empty_rows = "T";
    } else {
        $rm_empty_rows = "F";
    }

    # fetch row/col names
    my ($row_names, $col_names);
    my @row_names = map { $matrix->[$_]->[0] } (1 .. $#$matrix);
    my @col_names = @{$matrix->[0]}; shift @col_names;

    # cell 
    my $cnt = 0;
    my %scale = (0 => 1, 1 => 0.7, 2 => 0.4);
    my @bar_colors = ();
    my $bar_colors = "";
    my $alter_fun = <<ALTER;
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, (w-unit(0.15, "mm")) * $scale{0}, (h-unit(0.15, "mm")) * $scale{0}, 
        #grid.rect(x, y, w * $scale{0}, h * $scale{0},
            gp = gpar(fill = "$conf->{color_bg}", col = "$conf->{color_border}"))
    },
ALTER
    for my $gid (sort {$group->{$a}->{id} <=> $group->{$b}->{id}} keys %$group) {
        if($conf->{group_type} eq 'string') {

        $alter_fun .= <<ALTER;
    "$group->{$gid}->{name}" = function(x, y, w, h) {
        grid.rect(x, y, (w-unit(0.15, "mm")) * $scale{0}, (h-unit(0.15, "mm")) * $scale{$cnt}, 
        #grid.rect(x, y, w * $scale{0}, h * $scale{$cnt},
            gp = gpar(fill = "$conf->{color_group}->{$group->{$gid}->{name}}", 
            col = "$conf->{color_border}"))
    },
ALTER
        }
        else {

        $alter_fun .= <<ALTER;
    "$group->{$gid}->{name}" = function(x, y, w, h) {
        grid.rect(x, y, (w-unit(0.15, "mm")) * $scale{0}, (h-unit(0.15, "mm")) * $scale{0}, 
            gp = gpar(fill = "$conf->{color_group}->{$group->{$gid}->{name}}", 
            col = "$conf->{color_border}"))
    },
ALTER
        }
        push @bar_colors, "'$group->{$gid}->{name}' = ".
        "'$conf->{color_group}->{$group->{$gid}->{name}}'";

        $cnt ++;
        $cnt = 0 if($cnt > 2);
    }
    $alter_fun =~ s/,\n$/\n\)\n/;
    $bar_colors = "c(" . join(", ", @bar_colors) . ")";

    # bar
    my ($show_row_bar, $show_col_bar);
    if($conf->{show_row_bar} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        $show_row_bar = "rowAnnotation(rbar = anno_oncoprint_barplot())";
    } else {
        $show_row_bar = "NULL";
    }

    if($conf->{show_col_bar} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        $show_col_bar = "HeatmapAnnotation(cbar = anno_oncoprint_barplot())";
    } else {
        $show_col_bar = "NULL";
    }

    # title
    my ($row_title, $col_title, $legend_title);
    if($conf->{row_title} ne "") {
        $row_title = <<TMP;
row_title = "$conf->{row_title}",
    row_title_side = "$conf->{row_title_side}",
    row_title_gp = gpar(fontsize=$conf->{row_title_size}, fontface="$conf->{row_title_font}"),
TMP
    } else {
        $row_title = "";
    }

    if($conf->{col_title} ne "") {
        $col_title = <<TMP;
column_title = "$conf->{col_title}",
    column_title_side = "$conf->{col_title_side}",
    column_title_gp = gpar(fontsize=$conf->{col_title_size}, fontface="$conf->{col_title_font}"),
TMP
    } else {
        $col_title = "";
    }

    $legend_title = <<TMP;
heatmap_legend_param = list(
        title = "$conf->{legend_title}", 
        title_gp = gpar(fontsize=$conf->{legend_title_size}, fontface="$conf->{legend_title_font}"), 
        title_position = "$conf->{legend_title_pos}",
        legend_direction = "$conf->{legend_direction}"
    ),
TMP

    # row names | col names 
    my ($show_pct, $show_row_names, $show_col_names);
    if($conf->{show_pct} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        $show_pct = <<TMP;
show_pct = T,
    pct_gp = gpar(fontsize = $conf->{pct_size}, fontface = "$conf->{pct_font}"),
    pct_side = "$conf->{pct_side}",
    pct_digits = $conf->{pct_digits},
TMP
    } else {
        $show_pct = <<TMP;
show_pct = F,
TMP
    }

    if($conf->{show_row_names} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        $show_row_names = <<TMP;
show_row_names = T,
    row_names_gp = gpar(fontsize = $conf->{row_names_size}, fontface = "$conf->{row_names_font}"),
    row_names_side = "$conf->{row_names_side}",
    row_names_rot = $conf->{row_names_rot},
TMP
    } else {
        $show_row_names = <<TMP;
show_row_names = F,
TMP
    }

    if($conf->{show_col_names} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        $show_col_names = <<TMP;
show_column_names = T,
    column_names_gp = gpar(fontsize = $conf->{col_names_size}, fontface = "$conf->{col_names_font}"),
    column_names_side = "$conf->{col_names_side}",
    column_names_rot = $conf->{col_names_rot},
TMP
    } else {
        $show_col_names = <<TMP;
show_column_names = F,
TMP
    }

    my $rsc = <<RSC;
library(ComplexHeatmap)

mat_file <- "$matrixF"
width = $conf->{width}
height = $conf->{height}
row_order = as.logical("$conf->{keep_row_order}")
col_order = as.logical("$conf->{keep_col_order}")

mat <- read.table(mat_file, sep = "\\t", header = T, row.names = 1, check.names = F)
mat <- as.matrix(mat)
mat[is.na(mat)] <- ""
mat.rown <- rownames(mat)
mat.coln <- colnames(mat)

if(row_order) {
    row_order = mat.rown
} else {
    row_order = NULL
}
if(col_order) {
    col_order = mat.coln
} else {
    col_order = NULL
}

$alter_fun
col = $bar_colors

hm <- oncoPrint(
    mat, 
    get_type = function(x) strsplit(x, ";")[[1]],
    alter_fun = alter_fun,
    col = col,

    top_annotation = $show_col_bar,
    right_annotation = $show_row_bar,

    $row_title
    $col_title
    $legend_title

    $show_pct
    $show_row_names
    $show_col_names

    row_order = row_order, 
    column_order = col_order,
    remove_empty_columns = $rm_empty_cols,
    remove_empty_rows = $rm_empty_rows
)

pdf("$outpre.heatmap.pdf", width = width, height = height)
draw(hm, heatmap_legend_side = "$conf->{legend_side}")
dev.off()

RSC

    open my $ofh, "> $outpre.heatmap.R" or Error('can not write file', "$outpre.heatmap.R");
    print $ofh $rsc;
    close $ofh;

    chdir($outdir);
    system("$Rscript $outpre.heatmap.R");
    system("$convert -density 300 $outpre.heatmap.pdf $outpre.heatmap.png");
    chdir;

}

sub init_conf {

    my ($self, %opts) = @_;
    for my $opt (keys %opts) {
        if($opts{$opt}) {
            $self->{$opt} = $opts{$opt};
        }
    }

    return $self;
}

sub check_conf {

    my $self = shift;
    
    if(!-s $self->{matrix}) {
        Error('check_conf', "no matrix file");
    }
    if(!-s $self->{group}) {
        Error('check_conf', "no group file");
    }
    if(!defined $self->{outpre}) {
        Error('check_conf', "no outprefix setting");
    }
    else {
        $self->{outpre} = rel2abs($self->{outpre});
    }

    my $matrix = read_matrix($self->{matrix}, $conf, "F");
    my $group = read_group($self->{group});
    my @colors = split /[,;]/, $self->{color_list};

    if(scalar(@colors) < scalar(keys %$group)) {
        Error('check_conf', 'the colors cannot match the group');
    }
    for my $i (sort {$group->{$a}->{id} <=> $group->{$b}->{id}} keys %$group) {
        my $color = $group->{$i}->{color} ? $group->{$i}->{color} : 
                    $colors[$group->{$i}->{id}];
        $self->{color_group}->{$group->{$i}->{name}} = $color;
    }
    if($self->{color_border} eq 'NA' or $self->{color_border} eq '') {
        $self->{color_border} = "NA";
    }

    # inherit options
    for my $opt (qw/rm_empty_cols rm_empty_rows/) {
        if($self->{$opt} =~ /default/) {
            $self->{$opt} = $self->{rm_empty};
        }
    }
    for my $opt (qw/show_row_bar show_col_bar/) {
        if($self->{$opt} =~ /default/) {
            $self->{$opt} = $self->{show_bar};
        }
    }
    for my $opt (qw/row_title_size col_title_size legend_title_size/) {
        if($self->{$opt} =~ /default/) {
            $self->{$opt} = $self->{title_size};
        }
    }
    for my $opt (qw/pct_size row_names_size col_names_size/) {
        if($self->{$opt} =~ /default/) {
            $self->{$opt} = $self->{name_size};
        }
    }

    my $mat_col_num = scalar(@{$matrix->[0]}) - 1;
    my $mat_row_num = scalar(@$matrix) - 1;
    if($conf->{width} =~ /^auto$|^a$/i) {
        if($mat_col_num <= 10) {
            $conf->{width} = 6;
        } elsif($mat_col_num <= 20) {
            $conf->{width} = 8;
        } elsif($mat_col_num <= 30) {
            $conf->{width} = 10;
        } elsif($mat_col_num <= 100) {
            $conf->{width} = 10 + 8 * ($mat_col_num - 30) / 70;
        } else {
            $conf->{width} = 20;
        }
    }
    if($conf->{height} =~ /^auto$|^a$/i) {
        if($mat_row_num <= 10) {
            $conf->{height} = 6;
        } elsif($mat_row_num <= 20) {
            $conf->{height} = 8;
        } elsif($mat_row_num <= 30) {
            $conf->{height} = 10;
        } elsif($mat_row_num <= 100) {
            $conf->{height} = 10 + 8 * ($mat_row_num - 30) / 70;
        } else {
            $conf->{height} = 20;
        }
    }

    return $self;
}

sub normalize_matrix {

    my $matrix = shift;
    my $group  = shift;
    my $outF   = shift;
    my $type   = shift || 'auto';      # auto(a) | string(s) | numeric(n)

    my $ftype  = "NA";
    if($type eq 'auto' or $type eq 'a') {
        my $flag = 0;
        my ($p1, $p2, $m1, $m2) = (0, 0, 0, 0);
        for my $p (sort keys %$group) {
            $p1 ++;
            if($p =~ /[<=>]/) {
                $p2 ++;
            }
        }
        $flag += 1 if($p1 == $p2);

        for my $i (1 .. (scalar(@{$matrix->[0]}) - 1)) {
            for my $j (1 .. $#$matrix) {
                $matrix->[$j]->[$i] =~ /^\s+$/ and $matrix->[$j]->[$i] = '';
                $flag == 0 and $matrix->[$j]->[$i] =~ /^0$/ and $matrix->[$j]->[$i] = '';
                if($matrix->[$j]->[$i] ne 'NA' or $matrix->[$j]->[$i] ne '') {
                    $m1 ++;
                    if($matrix->[$j]->[$i] =~ 
                    /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
                        $m2 ++;
                    } else {
                        if($flag == 1) {
                            Warn("the matrix has unknown field: [$j, $i] => [$matrix->[$j]->[$i]]");
                        }
                    }
                }
            }
        }
        $flag += 2 if($m1 == $m2);

        if($flag == 0) {
            $ftype = 'string';
            ## donot remember the function of this code, but it dirive to another bug
            #for my $i (1 .. (scalar(@{$matrix->[0]}) - 1)) {
            #    for my $j (1 .. $#$matrix) {
            #        !defined $group->{$matrix->[$j]->[$i]} and $matrix->[$j]->[$i] = '';
            #    }
            #}
        }
        elsif($flag == 1) {
            Warn("the matrix format type is 'string' and the group format type is 'numeric'");
            $ftype = 'none';
        }
        elsif($flag == 2) {
            Warn("the matrix format type is 'numeric' and the group format type is 'string'");
            $ftype = 'none';
        }
        elsif($flag == 3) {
            $ftype = 'numeric';
        }
    }
    elsif($type eq 'string' or $type eq 's') {
        $ftype = 'string';
    }
    elsif($type eq 'numeric' or $type eq 'n') {
        $ftype = 'numeric';
    }

    my @fmatrix = ();
    if($ftype eq 'none') {
        Error("normalize_matrix", "the format types are not the same and the script would be exit");
        exit 1;
    }
    elsif($ftype eq 'string') {
        for my $i (0 .. (scalar(@{$matrix->[0]}) - 1)) {
            for my $j (0 .. $#$matrix) {
                if($i == 0 or $j == 0) {
                    $fmatrix[$j][$i] = $matrix->[$j]->[$i];
                }
                elsif($matrix->[$j]->[$i] ne 'NA' or $matrix->[$j]->[$i] ne '') {
                    my @seg = split /\s*[;,|]\s*/, $matrix->[$j]->[$i];
                    @seg = map { "'$group->{$_}->{name}'" } @seg;
                    $fmatrix[$j][$i] = join(";", @seg);
                }
                else {
                    $fmatrix[$j][$i] = "";
                }
            }
        }
    }
    elsif($ftype eq 'numeric') {
        for my $i (0 .. (scalar(@{$matrix->[0]}) - 1)) {
            for my $j (0 .. $#$matrix) {
                if($i == 0 or $j == 0) {
                    $fmatrix[$j][$i] = $matrix->[$j]->[$i];
                }
                elsif($matrix->[$j]->[$i] ne 'NA' or $matrix->[$j]->[$i] ne '') {
                    my $keys = expression_parse($matrix->[$j]->[$i], [sort keys %$group]);
                    my @class = map { $_ eq 'NA' ? 'NA' : "'$group->{$_}->{name}'" } @$keys;
                    my $class = join(";", @class);
                    if($class =~ /\bNA\b/) {
                        $fmatrix[$j][$i] = "";
                    }
                    else {
                        $fmatrix[$j][$i] = $class;
                    }
                }
                else {
                    $fmatrix[$j][$i] = "";
                }
            }
        }
    }

    my $len = 0;
    open my $ofh, "> $outF" or Error('cannot write file', $outF);
    for my $row (@fmatrix) {
        my $l = scalar(@$row);
        if($len == 0) {
            $len = $l;
        }
        elsif($len != $l) {
            Error('normalize_matrix', 'the number of columns is not equal');
        }
        print $ofh join("\t", @$row), "\n";
    }
    close $ofh;

    return $ftype;
}

sub read_matrix {

    my $matrixF = shift;
    my $conf = shift;
    my $verbose = shift || "T";

    my @matrix;
    my ($row_num, $col_num) = (0, 0);
    open my $fh, "< $matrixF" or Error("no matrix file", $matrixF);
    while (<$fh>) {
        chomp;
        $_ =~ s/\r$//;

        my @arr = ();
        if($. == 1) {
            @arr = split /\t/, $_;
        }
        else {
            @arr = split /\t/, $_, $col_num;
        }

        my $len = scalar(@arr);
        if($col_num == 0) {
            $col_num = $len;
        }
        elsif($col_num != $len) {
            Error('read_matrix', "the number of columns is not equal at line $. ".
            "<$col_num <=> $len>");
        }

        $row_num ++;
        if($conf->{limit} =~ /^\d+$/) {
            if($row_num <= $conf->{limit} + 1) {
                push @matrix, [@arr];
            }
        } else {
            push @matrix, [@arr];
        }
    }
    close $fh;
    if($verbose eq "T" and $conf->{limit} =~ /^\d+$/ and $row_num > $conf->{limit} + 1) {
        Warn("The row number of the matrix is [$row_num] and the limit is [$conf->{limit}]");
    }

    if(defined $conf->{transpose} and $conf->{transpose} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        my $tmatrix = transpose(\@matrix);
        @matrix = @$tmatrix;
    }

    return \@matrix;
}

sub read_group {
    my $groupF = shift;
    
    my %group = ();
    my $cnt = 0;
    open my $fh, "< $groupF" or Error("no group file", $groupF);
    while (<$fh>) {
        chomp;
        $_ =~ s/\r$//;
        $. == 1 and next;   # omit header
        /^#/ and next;
        /^$/ and next;
        my @arr = split /\t/, $_;
        if(scalar(@arr) == 2) {
            $group{$arr[0]}{id}    = $cnt;
            $group{$arr[0]}{name}  = $arr[1];
            $group{$arr[0]}{color} = "";
            $cnt ++;
        }
        elsif(scalar(@arr) == 3) {
            $group{$arr[0]}{id}    = $cnt;
            $group{$arr[0]}{name}  = $arr[1];
            $group{$arr[0]}{color} = $arr[2];
            $cnt ++;
        }
        else {
            Error("read_group", "bad format found at line $. .");
        }
    }
    close $fh;
    
    return \%group;
}

sub expression_parse {
    # we assume that the factor is 'x', and the compare symbols are: 
    # <=, <, =, >, >=
    my $x   = shift;
    my $exp = shift;
    
    my @class = ();
    for my $e (@$exp) {
        my $tmp = $e; $tmp =~ s/\s+//g;
        my @tmp = split /n|x/, $tmp;
        my @assert = ();
        for my $a (@tmp) {
            push @assert, $a if($a ne "")
        }
        if(scalar(@assert) == 0 or scalar(@assert) > 2) {
            Error("expression_class", "'$e' is not allowed");
        }

        my $flag = 0;
        for my $a (@assert) {
            if($a =~ /([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?/) {
                if(($` eq '<=' or $` eq '≤') and $x <= $&) {
                    $flag ++;
                }
                elsif(($` eq '<' or $` eq '＜') and $x < $&) {
                    $flag ++;
                }
                elsif(($` eq '=' or $` eq '==' or $` eq '＝') and $x == $&) {
                    $flag ++;
                }
                elsif(($` eq '>' or $` eq '＞') and $x > $&) {
                    $flag ++;
                }
                elsif(($` eq '>=' or $` eq '≥') and $x >= $&) {
                    $flag ++;
                }
                elsif(($' eq '<=' or $' eq '≤') and $& <= $x) {
                    $flag ++;
                }
                elsif(($' eq '<' or $' eq '＜') and $& < $x) {
                    $flag ++;
                }
                elsif(($' eq '=' or $' eq '==' or $' eq '＝') and $& == $x) {
                    $flag ++;
                }
                elsif(($' eq '>' or $' eq '＞') and $& > $x) {
                    $flag ++;
                }
                elsif(($' eq '>=' or $' eq '≥') and $& >= $x) {
                    $flag ++;
                }
            }
        }
        if($flag == scalar(@assert)) {
            push @class, $e;
        }
    }
    
    if(scalar(@class) == 0) {
        push @class, 'NA';
    }

    return \@class;
}

sub load_conf {

    my $confF = shift;
    my %opts = @_;

    if (!defined $confF or $confF eq "" or !-s $confF) {
        Error("load_conf", "no conf path: $confF");
    }
    
    timeLOG("found conf file: $confF");
    my $conf = LoadFile($confF);
    return $conf;
}

sub timeLOG {

    my $timeH = "[" . ftime() . "]";
    print "$timeH @_\n";
}

sub ftime {

    my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst) = localtime();
    $year += 1900;
    $mon += 1;

    my $ftime = sprintf("%d\-%02d\-%02d %02d:%02d:%02d", 
    $year,$mon,$day,$hour,$min,$sec);
    return $ftime;
}

sub Warn {

    print STDERR "[Warn]: @_\n";
}

sub Error {

    my $error = shift;
    my $ename = $error;
    print STDERR "[Error]: $error, [@_]\n";
    exit 1;
}

sub preset_color {

    my $num = shift;

    my @preset = ("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", 
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85");

    if($num <= 0 or $num > scalar(@preset)) {
        Error("preset_color", "the num [$num] is out of bounds");
    }
    my @set = @preset[0 .. ($num-1)];
    return [@set];
}

sub transpose {

    my $matrix = shift;

    my @matrix = @$matrix;
    my @tmatrix = map {
        my $x = $_;
        [ map { $matrix[$_][$x] } 0 .. $#matrix ];
    } 0 .. $#{$matrix[0]}; 

    return \@tmatrix;
}

sub clean_tmp_file {

    my $conf = shift;
    if($conf->{clean} =~ /^T$|^True$|^Y$|^Yes$|^1$/i) {
        my $outdir  = dirname($conf->{outpre});
        unlink("$outdir/Rplots.pdf") if(-s "$outdir/Rplots.pdf");
        unlink("$conf->{outpre}.heatmap.R") if(-s "$conf->{outpre}.heatmap.R");
        unlink("$conf->{outpre}.fmt.matrix.txt") if(-s "$conf->{outpre}.fmt.matrix.txt");
    }

    return 0;
}

sub usage {

    print <<HELP;

Usage:

set para in conf: 
    perl $0 --conf oncoheatmap.yml

set para in cmd:  
    perl $0 -m matrix.txt -g group.txt -o outprefix -l "red,blue"

Options:

HELP

}

