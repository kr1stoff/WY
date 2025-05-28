package GeneralPlot::Data;
#--------------------------------------------------#
#    Author:          Liu yubin                    |
#    Created time:    2020.05.07 15:48:05          |
#--------------------------------------------------#
use strict;
use warnings;
no strict 'refs';
use File::Basename;
use FindBin qw($Bin $RealBin);
use File::Spec::Functions qw(rel2abs);
use Config::General;
use parent 'Exporter';
our @EXPORT = qw(data_check data_deal);
our @EXPORT_OK = qw();

use lib "$Bin";
use lib "$Bin/..";

use GeneralPlot::Debug;

#-------------------------------------------------------------------------------
#  check data format and count needed dims
#-------------------------------------------------------------------------------
sub data_check {

    my $file = shift;
    my $type = shift;
    my $subfix = shift;
    my $sub = "format_$type";

    my @return = $sub->($file, $subfix);
    if(scalar(@return) == 0) {
        ERROR("data_check", "no return in data_check(), check the input file. ");
    }
    elsif(scalar(@return) == 1 and $return[0] == 0) {
        ERROR("data_check", "no data return in data_check(), check the input file. ");
    }

    return @return;
}

#-------------------------------------------------------------------------------
#  check data format and count needed dims
#-------------------------------------------------------------------------------
sub data_deal {
    
    my $file = shift;
    my $out  = shift;
    my $type = shift;
    
    my $sub = "deal_with_$type";
    my @return = $sub->($file, $out);    

}

#-------------------------------------------------------------------------------
#  format table
#-------------------------------------------------------------------------------
# table is like:
# id    attr1
# sample1   v1
# sample1   v2
# sample2   v3
# sample2   v4
#-------------------------------------------------------------------------------
sub format_table {

    my $file = shift;

    my %label1 = ();
    my %label2 = ();
    open my $fh, "< $file" or ERROR("format_table", "$!");
    while (<$fh>) {
        chomp;
        $. == 1 and next;
        my @arr = split /\t/, $_;
        $label1{$arr[0]} ++;
    }
    close $fh;

    my $label1_cnt = scalar(keys %label1);
    if($label1_cnt > 0) {
        return $label1_cnt;
    }
    else {
        ERROR("format_table", "data format wrong. ");
    }
}

#-------------------------------------------------------------------------------
#  format table2, call the first two columns table
#-------------------------------------------------------------------------------
# table2 is like:
# SNP   CHR POS P
# rs1   1   1   0.9148060
# rs2   1   2   0.9370754
# rs3   1   3   0.2861395
#------------------------------------------------------------------------------
sub format_table2 {

    my $file = shift;

    my %label1 = ();
    my %label2 = ();
    open my $fh, "< $file" or ERROR("format_table2", "$!");
    while (<$fh>) {
        chomp;
        $. == 1 and next;
        my @arr = split /\t/, $_;
        $label1{$arr[0]} ++;
        $label2{$arr[1]} ++;
    }
    close $fh;

    my $label1_cnt = scalar(keys %label1);
    my $label2_cnt = scalar(keys %label2);
    if($label1_cnt > 0 and $label2_cnt > 0) {
        #print STDERR "checking data ... \n";
        #print STDERR "column1: $label1_cnt; column2: $label2_cnt \n";
        return($label1_cnt, $label2_cnt);
    }
    else {
        ERROR("format_table2", "data format wrong. ");
    }
}

#-------------------------------------------------------------------------------
#  format table3, two or three columns table
#-------------------------------------------------------------------------------
# table3 is like:
# id    attr1   type
# sample1   v1  t1
# sample1   v2  t1
# sample2   v3  t2
# sample2   v4  t2
#-------------------------------------------------------------------------------
sub format_table3 {

    my $file = shift;

    my %label1 = ();
    my %label2 = ();
    open my $fh, "< $file" or ERROR("format_table3", "$!");
    while (<$fh>) {
        chomp;
        $. == 1 and next;
        my @arr = split /\t/, $_;
        $label1{$arr[0]} ++;

        if(scalar(@arr) == 3) {
            $label2{$arr[2]} ++;
        }
    }
    close $fh;

    my $label1_cnt = scalar(keys %label1);
    my $label2_cnt = scalar(keys %label2);
    if($label1_cnt > 0 and $label2_cnt > 0) {
        return($label1_cnt, $label2_cnt);
    }
    elsif($label1_cnt > 0) {
        return $label1_cnt;
    }
    else {
        ERROR("format_table3", "data format wrong. ");
    }
}

#-------------------------------------------------------------------------------
#  format matrix
#-------------------------------------------------------------------------------
# matrix is like:
# id    sample1    sample2    sample3
# attr1 v1  v2  v3
# attr2 v4  v5  v6
#-------------------------------------------------------------------------------
sub format_matrix {

    my $file = shift;
    my $subfix = shift || "none";

    my ($nrow, $ncol, $last_ncol) = (0, 0, 0);
    if($subfix eq 'none') {
        open my $fh, "< $file" or ERROR("format_matrix", "$!");
        while (<$fh>) {
            chomp;
            /^#/ and next;
            my @arr = split /\t/, $_;
            $nrow ++;
            $ncol = scalar(@arr);

            if($last_ncol == 0) {
                $last_ncol = $ncol;        
            }
            elsif($last_ncol != $ncol) {
                ERROR("format_matrix", "the number of column isnot equal at [$.]. ");
            }
        }
        close $fh;
    }
    else {
        open my $fh, "< $file" or ERROR("format_matrix", "$!");
        while (<$fh>) {
            chomp;
            my @arr = split /\t/, $_;
            if($. == 1) {
                for my $i (0 .. $#arr) {
                    if($arr[$i] =~ /$subfix$/) {
                        $ncol ++;
                    }
                }
            }
            $nrow ++;
        }
        close $fh;
    }

    return ($nrow, $ncol); 

}

#-------------------------------------------------------------------------------
#  format list
#-------------------------------------------------------------------------------
# list is like: 
# Group1    Group2
# v1    v1
# v2    v3
# v5    v6
# v7    v8
# the format list data is regular
#-------------------------------------------------------------------------------
sub format_list {

    my $file = shift;

    my ($nrow, $ncol, $last_ncol) = (0, 0, 0);

    open my $fh, "< $file" or ERROR("format_list", "$!");
    while (<$fh>) {
        chomp;
        /^#/ and next;
        my @arr = split /\t/, $_;
        $nrow ++;
        $ncol = scalar(@arr);

        if($last_ncol == 0) {
            $last_ncol = $ncol;
        
        }
        elsif($last_ncol != $ncol) {
            ERROR("format_list", "the number of column isnot equal at [$.]. ");
        }
    }
    close $fh;

    return ($ncol); 

}

#-------------------------------------------------------------------------------
#  format list2, for sbv venn 
#-------------------------------------------------------------------------------
# list2 is like :
# Group1    v1,v2,v3
# Group2    v1,v3,v4,v5
#-------------------------------------------------------------------------------
sub format_list2 {

    my $file = shift;
    
    my $count = 0;
    open my $fh, "< $file" or ERROR("format_list2", "$!");
    while (<$fh>) {
        chomp;
        /^#/ and next;
        /^$/ and next;
        $count ++;
    }
    close $fh;

    return $count;

}

#-------------------------------------------------------------------------------
#  format list3, for sbv venn
#-------------------------------------------------------------------------------
# list3 is like :
# Group1    Group2
# v1    v1
# v2    v3
# v5    v6
#       v8
#-------------------------------------------------------------------------------
sub format_list3 {

    my $file = shift;

    my $count = 0;
    my $line_cnt = 0;
    open my $fh, "< $file" or ERROR("format_list3", "$!");
    while (<$fh>) {
        chomp;
        /^#/ and next;
        /^$/ and next;
        if($line_cnt == 0) {
            my @arr = split /\t/, $_;
            $count = scalar(@arr);
            $line_cnt = 1;
        }
    }
    close $fh;

    return $count;

}

#-------------------------------------------------------------------------------
#  format group, call the first two columns in no header table
#-------------------------------------------------------------------------------
# group is like:
# id1   group1
# id2   group1
# ...
#------------------------------------------------------------------------------
sub format_group {

    my $file = shift;

    my %label1 = ();
    my %label2 = ();
    open my $fh, "< $file" or ERROR("format_group", "$!");
    while (<$fh>) {
        chomp;
        my @arr = split /\t/, $_;
        $label1{$arr[0]} ++;
        $label2{$arr[1]} ++;
    }
    close $fh;

    my $label1_cnt = scalar(keys %label1);
    my $label2_cnt = scalar(keys %label2);
    if($label1_cnt > 0 and $label2_cnt > 0) {
        return($label1_cnt, $label2_cnt);
    }
    else {
        ERROR("format_group", "data format wrong. ");
    }
}


#-------------------------------------------------------------------------------
#  deal with the data in format list3, fill with "NA"
#-------------------------------------------------------------------------------
sub deal_with_list3 {

    my $file = shift;
    my $out  = shift;

    open my $fh, "< $file" or ERROR("deal_with_list3", "$!");
    open my $ofh, "> $out" or ERROR("deal_with_list3", "$!");
    while (<$fh>) {
        chomp;
        my @arr = split /\t/, $_;
        if($. == 1) {
            for my $v (@arr) {
                if($v eq '') {
                    ERROR("deal_with_list3", "sample id '' is not regular. ");
                }
                elsif($v eq 'NA') {
                    ERROR("deal_with_list3", "sample id 'NA' is not regular. ");
                }
            }
        }
        @arr = map { $_ eq "" ? "NA" : $_ } @arr;
        print $ofh join("\t", @arr), "\n";
    }
    close $fh;
    close $ofh;
    
    return 0;

}

#-------------------------------------------------------------------------------
#  format sankey
#-------------------------------------------------------------------------------
# sankey is like:
# attr1    attr2    attrn...    Freq
# A        C        ...         1
# A        D        ...         1
# B        C        ...         1
# the last column should be named as Freq
#-------------------------------------------------------------------------------
sub format_sankey {

    my $file = shift;

    timeLOG("format check for 'sankey' ...");
    my $head_len = 0;
    open my $fh, "< $file" or die "$!";
    while (<$fh>) {
        chomp;
        $_ =~ s/\r$//;
        my @arr = split /\t/, $_;
        if($. == 1) {
            $head_len = scalar(@arr);
            if($head_len < 3) {
                _error('format_sankey', "the number of columns should be equal or greater than 3");
            }
            if($arr[-1] ne 'Freq') {
                _error('format_sankey', "the last column should be named as 'Freq'");
            }
        }
        else {
            my $tmp_len = scalar(@arr);
            if($head_len != $tmp_len) {
                _error('format_sankey', "there are $head_len columns in header,", 
                "but $tmp_len in line $.");
            }
            if($arr[-1] !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
                _error('format_sankey', "the last column 'Freq' should be numeric,",
                "line $., column $head_len => '$arr[-1]'");
            }
        }
    }
    close $fh;
    timeLOG("format check ok.");

    return $head_len;
}

#-------------------------------------------------------------------------------
#  format polarbar
#-------------------------------------------------------------------------------
# the file is like:
# id   number
# A    1
# B    3
# C    7
# only two columns 
#-------------------------------------------------------------------------------
sub format_polarbar {

    my $file = shift;

    timeLOG("format check for 'polarbar' ...");
    my $head_len = 0;
    open my $fh, "< $file" or die "$!";
    while (<$fh>) {
        chomp;
        $_ =~ s/\r$//;
        my @arr = split /\t/, $_;
        if($. == 1) {
            $head_len = scalar(@arr);
            if($head_len != 2) {
                _error('format_polarbar', "the number of columns should be 2, header is $head_len");
            }
        }
        else {
            my $tmp_len = scalar(@arr);
            if($head_len != $tmp_len) {
                _error('format_polarbar', "there are $head_len columns in header,", 
                "but $tmp_len in line $.");
            }
            if($arr[-1] !~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
                _error('format_polarbar', "the last column should be numeric,",
                "line $., column $head_len => '$arr[-1]'");
            }
        }
    }
    close $fh;
    timeLOG("format check ok.");

    return $head_len;

}

#-------------------------------------------------------------------------------
#  format splitviolin
#-------------------------------------------------------------------------------
#  be same to matrix
#-------------------------------------------------------------------------------
sub format_splitviolin {

    my $file = shift;
    my $subfix = shift || "none";

    timeLOG("format check for 'splitviolin' ...");
    timeLOG("file: '$file' ");
    my ($nrow, $ncol, $head_len) = (0, 0, 0);
    if($subfix eq 'none') {
        open my $fh, "< $file" or die "$!";
        while (<$fh>) {
            chomp;
            $_ =~ s/\r$//;
            my @arr = split /\t/, $_;
            $nrow ++;
            $ncol = scalar(@arr);

            if($. == 1) {
                $head_len = scalar(@arr);
            }
            else {
                if($head_len != $ncol) {
                    _error('format_splitviolin', "there are $head_len columns in header,",
                    "but $ncol in line $.");
                }
            }
        }
        close $fh;
    }
    else {
        open my $fh, "< $file" or die "$!";
        while (<$fh>) {
            chomp;
            $_ =~ s/\r$//;
            my @arr = split /\t/, $_;
            if($. == 1) {
                for my $i (0 .. $#arr) {
                    if($arr[$i] =~ /$subfix$/) {
                        $ncol ++;
                    }
                }
                $head_len = scalar(@arr);
            }
            else {
                my $tmp_len = scalar(@arr);
                if($head_len != $tmp_len) {
                    _error('format_splitviolin', "there are $head_len columns in header,",
                    "but $tmp_len in line $.");
                }
            }
            $nrow ++;
        }
        close $fh;
    }
    timeLOG("format check ok.");

    return ($nrow, $ncol); 

}

sub _error {
    my $ename = shift;
    my $error = $ename;
    print STDERR "[Error]: $error, [@_]\n";
    exit 1;
}

return 1;

