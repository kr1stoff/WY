package SBV::FontDefault;

#-------------------------------------------------------------------------------
# Date: 09/04/2013 04:04:19 PM 
#-------------------------------------------------------------------------------

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(textWidth textHeight);

#use PostScript::Font::TTtoType42;
use PostScript::FontMetrics;
use File::Basename qw(basename);

use Math::Round;
use FindBin;

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/..";
use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/../lib";
use SBV;
use SBV::DEBUG;

my $allfonts = {};

# load font files
# for read fast, now just support afm files
# afm files can be extract from ttf file by script 'ttf2afm.pl'
sub load_font
{
	my @path = @_;

	foreach my$path (@path)
	{
		my @afm_files = _glob_files($path,'afm');
		#my @ttf_files = _glob_files($path,'ttf');
		#my @otf_files = _glob_files($path,'otf');
		_load_afm_files(@afm_files);
		#_load_ttf_files(@ttf_files);
	}
	
	return 1;
}

# get the files who has the specific suffix
sub _glob_files
{
	my $path = shift;
	my $suffix = shift;
	my @files;

	opendir DIR,$path or die "$path $!";
	while(my $filename = readdir(DIR))
	{
		push @files , "$path/$filename" if ($filename =~ /\.$suffix$/i)
	}
	closedir DIR;

	return @files;
}

# load ttf files
sub _load_ttf_files
{
	my @files = @_;
	my @afm_files; 
	foreach my$file (@files)
	{
		my $font = PostScript::Font::TTtoType42->new($file);
		$font->afm_as_string;
		my $weight = _weight_to_num($font->_str(2)) || 0;
		my $name = $font->_str(1);
		$allfonts->{$name}->{$weight} = $file;	
	}
}

# load afm files
sub _load_afm_files
{
	my @files = @_;

	foreach my$file(@files)
	{
		my $info = new PostScript::FontMetrics ($file);
		my $weight = _weight_to_num($info->weight) || 0;
		my $name = lc $info->FamilyName;
		$allfonts->{$name}->{$weight} = $file;
	}

	return 1;
}

# turn postscript weight to num 
sub _weight_to_num
{	
	my $weight = lc shift;
	my %hash = ('regular' => 0,'normal'=>0,'bold'=>1,'italic'=>2,'oblique'=>2,
		'bold italic'=>3,'bold oblique'=>3);
	
	if (exists $hash{$weight}){return $hash{$weight} }
	else{ WARN('the font weight is error:',$weight) }
}

# fetch the font hash from SBV::allStyle
sub fetch_font
{
	my $class = shift;
	my $name = shift || "default";
	my $font = {};
	
	$name = "CLASS$name" if ($name !~ /^CLASS/ && $name !~ /^ID/ && $name ne "default");

	my $default = $SBV::allStyle->{text}->{default};

	if (! exists $SBV::allStyle->{text}->{$name})
	{
		$font = $default;	
	}
	else
	{
		my $selfDef = $SBV::allStyle->{text}->{$name};
		
		my @tag = ("font-size","font-family","font-weight","font-style","fill");
		
		map {
			$font->{$_} = defined $selfDef->{$_} ? $selfDef->{$_} : $default->{$_};
		} @tag;

		$font->{'font-size'} = trans_size($font->{'font-size'}) . 'px';
	}
	
	check_font($font);
	
	bless $font , $class;

	return $font;
}

# create new font
sub new
{
	my $class = shift;
	my $style = shift;
	
	my @fonts_path = ("/home/aipeng/work/develepment/SBV/fonts/afm/");
	load_font(@fonts_path);

	my $font = {};
	my $default = $SBV::allStyle->{text}->{default};

	return $class->fetch_font() unless (defined $style);

	if (ref $style eq "HASH")
	{
		my @tag = ("font-size","font-family","font-weight","font-style","fill");
		map { $font->{$_} = defined $style->{$_} ? $style->{$_} : $default->{$_} } @tag;
	}
	else
	{
		my @attrs = split /[;|]/ , $style;
		foreach my$attr(@attrs)
		{
			my ($name,$val) = split /:/ , $attr;
			if ($name eq "fill" || $name eq "color")
			{
				$val = SBV::Colors::fetch_color($val);
				$name = "fill";
			}
			elsif ($name !~ /font/)
			{
				$name = "font-$name";	
			}
			
			$font->{$name} = $val;
		}

		my @tag = ("font-size","font-family","font-weight","font-style","fill");
		map { $font->{$_} = defined $font->{$_} ? $font->{$_} : $default->{$_} } @tag;
	}
	
	$font->{'font-size'} = trans_size($font->{'font-size'}) . 'px';
	
	check_font($font);

	bless $font , $class;
	return $font;
}

# check the font
sub check_font
{
	my $font = shift;
	my $ffam = lc $font->{'font-family'};
	
	if (! exists $allfonts->{$ffam})
	{
		WARN('the font family is not supported now, will use the Arail instead',$ffam);
		$font->{'font-family'} = 'arial';
	}

	return 1;
}

# generate font style string
sub toStyle
{
	my $self = shift;
	my $str = "";

	foreach my $key (keys %$self)
	{
		$str .= "$key\:$self->{$key};";	
	}

	return $str;
}

# set the attr for font 
sub setAttr
{
	my $self = shift;
	my $style = shift;

	if (ref $style eq "HASH")
	{
		foreach my$name(keys %$style)
		{	$name = "font-$name" if ($name ne "fill" && $name !~ /font/);
			$self->{$name} = $style->{$name};
		}
	}
	else
	{
		my @attrs = split /;/ , $style;
		foreach my$attr(@attrs)
		{
			my ($name,$val) = split /:/ , $attr;
			$name = "font-$name" if ($name ne "fill" && $name !~ /font/);
			$self->{$name} = $val;
		}
	}

	return $self;
}

# trans the size unit to px
sub trans_size
{
	my $fsize = shift;
	my ($num,$unit) = $fsize =~ /^(\d+)([a-z]*)$/;
	
	if ($unit eq "" || $unit eq "px")
	{
		return $num;
	}
	elsif ($unit eq "pt")
	{
		$num *= 4/3;
		return $num;
	}
}

# fetch the text width
sub fetch_text_width
{
	my $self = shift;
	my $str = shift;
	my $hscale = shift || 1;

	return nearest 0.01 , textWidth($self,$str,$hscale);
}

# fetch max text width
sub fetch_max_text_width
{
	my $self = shift;
	my $text = shift;
	my $hscale = shift || 1;

	my $max = textWidth($self,$text->[0],$hscale);

	map { my$len=textWidth($self,$_,$hscale); $max = $len if ($len > $max) } @$text;

	return nearest 0.01 , $max;
}

# fetch the text height
sub fetch_text_height
{
	my $self = shift;
	my $vscale = shift || 1;

	return nearest 0.01 , textHeight($self->{'font-size'},$vscale);
}

# return the height of a text
# the default unit is pt
# px = pt * 0.75
sub textHeight{
	my ($fsize,$vscale) = @_;
	$fsize = $1 if ($fsize =~ /^(\d+)/);
	$fsize = trans_size($fsize);
	$vscale = 1 if(!$vscale); # vertical scale

	my $factor = 0.73;
	my $h = nearest 0.001 , ($fsize * $factor);
	return $h + 2;
}

# returns the width of a text string, add by fanwei
# Modified by aipeng
sub textWidth
{
	my $font = shift;
	my $str = shift;
	my $hscale = shift || 1; #horizone scale

	# replace the tab in string by 4 spaces
	my $tab = "\t";
	my $replace = " " x 4;
	$str =~ s/$tab/$replace/g;
	
	my $fsize  = $font->{'font-size'};
	$fsize = trans_size($fsize) * $hscale;

	my $ffam   = lc $font->{'font-family'};
	my $fwei   = $font->{'font-weight'};
	my $fstyle = $font->{'font-style'};
	
	my $num = _weight_to_num($fwei) + _weight_to_num($fstyle);
	$num = 0 if (not exists $allfonts->{$ffam}->{$num});
	my $info = new PostScript::FontMetrics ($allfonts->{$ffam}->{$num});	
	
	return $info->kstringwidth($str,$fsize); 
}
