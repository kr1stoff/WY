package SBV::IMAGE::HEATMAP;
#------------------------------------------------+
#    [APM] This moudle is generated by amp.pl    |
#    [APM] Created time: 2013-10-14 16:58:58     |
#------------------------------------------------+
=pod

=head1 Name

SBV::IMAGE::HEATMAP -- moudle to draw heatmap figure for matrix data

=head1 Synopsis

This module is not meant to be used directly

=head1 Feedback

Author: Peng Ai
Email:  aipeng0520@163.com

=head1 Version

Version history

=head2 v1.0

Date: 2013-10-14 16:58:58

=cut

use strict;
use warnings;
require Exporter;

use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/..";
use lib "$FindBin::RealBin/../lib";

use SBV::IMAGE::MERGE;
use SBV::Hcluster;
use SBV::STAT;
use SBV::DEBUG;

sub new 
{
	my ($class,$file,$conf)	= @_;
	my $heatmap = {};
	
	$heatmap->{file} = $file;
	$heatmap->{conf} = $conf;
	timeLOG("start to read data ... ");
	my $data = SBV::DATA::Frame->new($file,header=>$conf->{header},rownames=>$conf->{rownames});
	timeLOG("read data done ... ");
	$heatmap->{data} = $data;
	bless $heatmap , $class;
	return $heatmap;
}

sub plot
{
	my ($self,$parent) = @_;

	my $conf = $self->{conf};
	my $file = $self->{file};
	my $data = $self->{data};
	
	my @names = $data->names;
	my @rownames = $data->rownames;
	
	# draw the background for whole heatmap
	SBV::DRAW::background($conf,$parent);
	my $group = $parent->group(id=>"heatmap$SBV::idnum");
	$SBV::idnum ++;
	
	# init the figure layout and size  
	my ($horizontal,$vertical,$main_width,$main_height) = _init_layout($conf,$data);
	timeLOG("init figure size done ... ");

	# create the temporary folder
	my $path = "HEATMAP" . time();
	mkdir $path;
	
	my ($hcluster,$vcluster);

	# regenerate the file if vertical is tree 
	if ($conf->{vertical}->{order} eq "tree" || $conf->{horizontal}->{order} eq "tree")
	{
		my $tempfile = $file;
		if (0 == $conf->{header})
		{
			$tempfile = SBV::DATA::add_header($file,"$path/temp.addheader",rownames=>$conf->{rownames});
		}
		my $cluster = SBV::Hcluster->new($tempfile);
		
		if ($conf->{horizontal}->{order} eq "tree")
		{
			timeLOG('start to do horizontal cluster ... ');
			$cluster->hcluster(transpose=>0);
			$cluster->save(file=>"$path/temp.horizontal.newick",name=>\@rownames);
			$hcluster = "$path/temp.horizontal.newick";
		}

		if ($conf->{vertical}->{order} eq "tree")
		{
			timeLOG('start to do vertical cluster ... ');
			$cluster->hcluster(transpose=>1);
			$cluster->save(file=>"$path/temp.vertical.newick",name=>\@names);
			$vcluster = "$path/temp.vertical.newick";
			
			my $treeio = Bio::TreeIO->new('-format'=>'newick',-file=>"$path/temp.vertical.newick");
			my $tree = $treeio->next_tree;
			my @leaves = map { $_->id } $tree->get_leaf_nodes;
			$data->save(file=>"$path/temp.dat",names=>\@leaves,header=>1);
			$file = "$path/temp.dat";
		}
	}
	
	# push the heatmap to horizontal as the first dataset
	my $heatset = {
		type => 'heatmap',
		file => $file,
		color => "$conf->{colors}",
		header => 1,
		width => $main_width,
		scale => $conf->{scale},
		header => $conf->{header},
		rownames => $conf->{rownames}
	};
	
	$heatset->{min} = $conf->{min} if (defined $conf->{min});
	$heatset->{max} = $conf->{max} if (defined $conf->{max});
	
	unshift @{$horizontal->{datasets}->{dataset}} , $heatset;
	
	# push the heatmap to vertical as the first dataset (no draw)
	my $pseudo_set = {
		type => 'heatmap',
		file => $file,
		color => "$conf->{colors}",
		header => 1,
		width => $main_height,
		sclae => $conf->{scale},
		show => 0
	};
	unshift @{$vertical->{datasets}->{dataset}} , $pseudo_set if ($vertical->{order} eq "tree");
	
	# draw horizontal map 
	if ($horizontal->{order} eq "tree")
	{
		$horizontal->{oriental} = "left";
		my $tree = SBV::IMAGE::TREE->new($hcluster,$horizontal);
		$tree->plot($group);
	}
	elsif ($horizontal->{order} eq "default")
	{
		my $merge = SBV::IMAGE::MERGE->new(\@rownames,$horizontal);
		$merge->plot($group,'horizontal');
	}
	elsif ($horizontal->{order} eq "taxonomy")
	{
		
	}
	else
	{
		ERROR('err_heatmap_order',$horizontal->{order});
	}
	
	# draw vertical map 
	if ($vertical->{order} eq "tree")
	{
		$vertical->{oriental} = "bottom";
		my $tree = SBV::IMAGE::TREE->new($vcluster,$vertical);
		$tree->plot($group);
	}
	elsif ($vertical->{order} eq "default")
	{
		my $merge = SBV::IMAGE::MERGE->new(\@names,$vertical);
		$merge->plot($group,'vertical');
	}
	else
	{
		ERROR('err_heatmap_order',$vertical->{order});
	}

	# remove the temporary file and path
	unless ($conf->{keep_temp_file})
	{
		unlink ("$path/temp.vertical.newick");
		unlink ("$path/temp.horizontal.newick");
		unlink ("$path/temp.dat");
		rmdir $path;
	}

	#--------------------------------------------------------------------
	# raw color bar legend for heatmap
	# init the axis of fpkm val
	my $hi = $SBV::conf->{hspace};
	
	my @max = map { max($data->{row}->{$_}) } @rownames;
	my $max = max(\@max);
	my @min = map { min($data->{row}->{$_}) } @rownames;
	my $min = min(\@min);
	$min = $conf->{min} if (defined $conf->{min});
	$max = $conf->{max} if (defined $conf->{max});
	my $tick = SBV::STAT::dividing($min,$max,-nture=>1,-xtrue=>1);
	$tick = "-1 1 1" if ($conf->{scale} ne "none");

	my $color_axis = SBV::STONE::AXIS->new(tick=>$tick,side=>"right",length=>$vertical->{tw},
		parent=>$group,bone=>0,translate=>0,start=>0,skip_first_tick=>0,skip_last_tick=>0);
	my $axis_thick = $color_axis->thickness;

	#my @colors = map { SBV::Colors::fetch_color($_) } split /[\s\t\,\;]+/ , $conf->{colors};
	my @colors = SBV::Colors::fetch_brewer_color($conf->{colors});
    my $fill = SBV::Colors::gradient(\@colors);

	my $ver_legend_size = 60;
	my $rectx = $vertical->{ox};
	my $recty = $conf->{oty} + 2*$hi;
	my $height = $ver_legend_size - 4*$hi - $axis_thick;
	$group->rect(x=>$rectx,y=>$recty,width=>$vertical->{tw},height=>$height,style=>"fill:$fill;stroke-width:0");
	$color_axis->plot(ox=>$rectx,oy=>$recty+$height);

	# draw legend
	if ($conf->{legend})
	{
		my $legend = SBV::STONE::LEGEND->new(conf=>$conf->{legend});
		$legend->location($conf);
		$legend->draw($parent);
	}
}

sub _init_layout
{
	my $conf = shift;
	my $data = shift;

	my $horizontal = $conf->{horizontal};
	my $vertical = $conf->{vertical};
	
	my $hi = $SBV::conf->{hspace};
	my $vi = $SBV::conf->{vspace};

	# the size of label and the tree( if defined)
	timeLOG("start to init the label size ... ");
	my $font = SBV::Font->fetch_font('leaf');
	my @names = $data->names;
	my @rownames = $data->rownames;
	
	my $hor_main_size = $conf->{cell_width} * ($#names+1);
	my $ver_main_size = $conf->{cell_height} * ($#rownames+1);
	
	# the size of datasets
	timeLOG("start to init the datasets size ... ");
	my $hor_datasets_size = _fetch_datasets_width($horizontal);
	my $ver_datasets_size = _fetch_datasets_width($vertical);
	
	# the size of the gradient color legend bar ( at the bottom of the figure height 
	my $ver_legend_size = 60;

	my $hor_label_width = 0;
	if ($horizontal->{hide_names})
	{
		$hor_label_width = 0;
	}
	elsif ($horizontal->{max_label_width})
	{
		$hor_label_width = $horizontal->{max_label_width};
	}
	else 
	{
		$hor_label_width = $font->fetch_max_text_width(\@rownames);
	}
	
	my $ver_label_width = $vertical->{hide_names} ? 0 : $font->fetch_max_text_width(\@names);
	
	my $hor_tree_size = $horizontal->{order} eq "tree" ? 
		$conf->{tw} - $hor_main_size - $hor_label_width - 4*$hi - $hor_datasets_size : 0;
	my $ver_tree_size = $vertical->{order} eq "tree" ? 
		$conf->{th} - $ver_main_size - $ver_label_width - 2*$vi - $ver_datasets_size - $ver_legend_size : 
		$conf->{th} - $ver_main_size - $ver_label_width - 2*$vi - $ver_datasets_size - $ver_legend_size ;
	
	my $hor_header_size = $hor_label_width + $hor_tree_size + 2*$hi;
	my $ver_header_size = $ver_label_width + $ver_tree_size + 2*$vi;

	# init the horizontal and vertical conf 
	$horizontal->{y} = $conf->{oy} - $ver_header_size;
	$horizontal->{height} = $ver_main_size;
	$vertical->{x} = $conf->{ox} + $hor_header_size;
	$vertical->{width} = $hor_main_size;
	$vertical->{y} = $conf->{oy};
	$vertical->{height} = $conf->{th} - $ver_legend_size;
	
	SBV::CONF::inherit($conf,$horizontal);
	SBV::CONF::inherit($conf,$vertical);
	
	# draw the split line 
	if ($conf->{show_split_line})
	{
		my $splitG = $SBV::svg->group();

		my $x0 = $conf->{ox};
		my $x1 = $x0 + $hor_header_size;
		my $x2 = $x1 + $hor_main_size;
		my $x3 = $conf->{ox} + $conf->{tw};
		my $y0 = $conf->{oty};
		my $y1 = $y0 + $ver_legend_size;
		my $y2 = $y1 + $ver_datasets_size;
		my $y3 = $y2 + $ver_main_size;
		my $y4 = $y3 + $ver_header_size;

		$splitG->line(x1=>$x0,x2=>$x3,y1=>$y1,y2=>$y1,class=>"split");
		$splitG->line(x1=>$x0,x2=>$x3,y1=>$y2,y2=>$y2,class=>"split");
		$splitG->line(x1=>$x0,x2=>$x3,y1=>$y3,y2=>$y3,class=>"split");

		$splitG->line(x1=>$x1,x2=>$x1,y1=>$y1,y2=>$y4,class=>"split");
		$splitG->line(x1=>$x2,x2=>$x2,y1=>$y1,y2=>$y4,class=>"split");
	}

	return ($horizontal,$vertical,$hor_main_size,$ver_main_size);
}

sub _fetch_datasets_width
{
	my $conf = shift;
	return 0 if (! $conf->{datasets});
	
	$conf = SBV::CONF::fetch_first_conf("datasets",$conf,0);
	return 0 if (! $conf->{dataset});

	my $width = 0;
	my $default_width = $conf->{width} || 16;
	my $vi = $SBV::conf->{vspace};
	if (ref $conf->{dataset} eq "ARRAY")
	{
		foreach my$subconf(@{$conf->{dataset}})
		{
			next if (! defined $subconf->{file});
			my $dataset_width = $subconf->{width} || $default_width;
			$width += $dataset_width + 2*$vi;
		}
	}
	elsif (ref $conf->{dataset} eq "HASH")
	{
		my $dataset_width = $conf->{dataset}->{width} || $default_width;
		$width += $dataset_width + 2*$vi;
	}
	$width += 2*$vi;

	return $width;
}
