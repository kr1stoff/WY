#!/Bio/bin/perl
#-----------------------------------------------+
#    [APM] This script was created by amp.pl    |
#    [APM] Created time: 2019-08-07 14:11:24    |
#-----------------------------------------------+
# name: auto_branch.pl
# func: 
# version: 1.0

use strict;
use warnings;

use Getopt::Std;
use Bio::TreeIO;
use POSIX;
use File::Basename qw/dirname basename/;
use Config::General;
use Data::Dumper;

my %opts = (r=>0.6,t=>"clade",i=>0,m=>"circular",b=>60,f=>"simple");
getopts('r:t:i:m:f:',\%opts);

die qq(
Usage:   perl $0 <tree_file> <leaves color def file>

Options: -r FLOAT    the ratio of leaves to color a clade, [0.6]
         -t STR      the type of color clade, 'clade' line or clade 'range', [clade]
         -m STR      the tree type, 'circular' or 'unrooted', [circular]
         -i INT      ignore the branch length, 1 means yes, 0 means not, [0]
         -b INT      bootstrap threshold for -i==1, [60]
         -f STR      the format of <leaves color def file>, 'simple' or 'detail'

Note: 
<leaves color def file> simple format (separated by <TAB>):
leaf_name\tcolor

<leaves color def file> detail format (separated by <TAB>):
order\tclient_id\tvcf_id\tgroup\tcolor

) unless @ARGV == 2;

my $tree_file  = shift @ARGV;
my $color_file = shift @ARGV;
my $ratio = $opts{r};
my $group_num = 0;

# fetch the name
my $name = basename($tree_file);
$name =~ s/\.(\w+)$//;

# read the tree file
my $treeio = Bio::TreeIO->new('-format'=>"newick",-file=>$tree_file);
my $tree = $treeio->next_tree;
my $root_node = $tree->get_root_node;

# read the leaves color info
my %colors = read_leves_color($color_file);

# auto calc the clade color with 'most' 
my %clade_colors;
calc_clade_color($root_node,\%clade_colors);

open my $ofh , ">clade.def" or die $!;
foreach (sort {$b<=>$a} keys %clade_colors){
    print $ofh "$_\t$opts{t}\t$clade_colors{$_}\n";
}

foreach ( keys %colors ){
    print $ofh "$_\t$opts{t}\t$colors{$_}\n";
}
close $ofh;

# set the configuration file for sbv
my $sbv = "perl /home/aipeng/work/develepment/SBV/bin/sbv.pl";
my $demo_conf = "/home/aipeng/work/develepment/SBV/demo/circular_tree/circular_tree.conf";
my $etc_path  = "/home/aipeng/work/develepment/SBV/etc";

my $config = Config::General->new(-ConfigFile=>$demo_conf,-ConfigPath=>[$etc_path],
    -LowerCaseNames=>1,-IncludeAgain=>1,-AutoTrue=>1,-SplitPolicy=>'equalsign');
my %config = $config->getall;

$config{tree}{model} = $opts{m};
$config{tree}{ignore_branch_length} = $opts{i};
$config{tree}{bootstrap}{text} = $opts{i};
$config{tree}{bootstrap}{threshold} = $opts{b} if ($opts{i});
$config{tree}{linkage_type} = "line" if ($opts{i});
$config{legends} = {} if 1 == $group_num;

$config->save_file("$name.tree.conf",\%config);

# print Dumper($config{styles}{text});

# run sbv to plot treee
system("$sbv tree -conf $name.tree.conf -out ${name}.$opts{m}_tree $tree_file");

#-------------------------------------------------------------------------------
#  sub functions
#-------------------------------------------------------------------------------
sub read_leves_color {
    my $file = shift;
    my %colors;

    open my $fh_color , $file or die $!;
    open my $ofh , ">legends.txt" or die $!;
    
    my %groups;

    <$fh_color> if ($opts{f} eq "detail");

    while(<$fh_color>){
        chomp;
        next unless $_;
        my ($leaf,$color,$group);
        if ($opts{f} eq "simple"){
            ($leaf,$color,$group) = (split /\t/)[0,1,2];
        }else{
            ($leaf,$group,$color) = (split /\t/)[1,3,4];
        }
        $color =~ s/^#//;
        $colors{$leaf} = $color;
        
        unless ($groups{$group}){
            my $legend_str = $opts{t} eq "clade" ? "$group\tshape=17;stroke_width=2;color=$color\n" : "$group\tshape=42;fill=$color\n";
            print $ofh $legend_str;
            $groups{$group} = 1;
            $group_num ++;
        }
    }
    close $fh_color;
    close $ofh;

    return %colors;
}

sub calc_clade_color {
    my $root   = shift;
    my $clade_colors = shift;
    return if $root->is_Leaf;
    
    my @nodes = $root->each_Descendent;
    foreach my $node (@nodes) {
        calc_clade_color($node,$clade_colors);
    }
    
    my @leaf_colors = map { my$id = $_->id; $colors{$id} or die "the color of [$id] is not defined" } grep { $_->is_Leaf } $root->get_all_Descendents;
    my $root_color  = _calc_clade_color(@leaf_colors);
    $clade_colors->{$root->internal_id} = $root_color if ($root_color);
}

sub _calc_clade_color {
    my @colors = @_;
    my $sum = scalar @colors;
    my $cutoff = ceil ( $sum * $ratio );

    my %counts;
    for ( @colors ){ $counts{$_} ++; }
    for ( keys %counts ) {  return $_ if $counts{$_} >= $cutoff } 

    return;
}
