use strict;

my $table = shift;
my $blast = shift;

my %hash;
my ($name,$md5);
open A,$table || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    ($name,$md5) = (split /\t/,$_)[0,1];
    $hash{$md5} = $name;
}
close A;

open B,$blast || die $!;
while(<B>){
    chomp;
    my @line = split /\t/,$_;
    my $md = $line[0];
    my $other = join "\t",@line[1..$#line];
    my $output = `wc -l $blast`;
    my ($line_count) = $output =~ /(\d+)/;
    if ($line_count == 1){
        if ($_ =~ /#/){
            print "$_\n";
        }
        print "$name\t0\t-\n";
    }
    else{
        if ($_ =~ /#/){
            print "$_\n";
        }
        else{
            if(exists $hash{$md}){
                print "$hash{$md}\t$other\n";
            }
            else{
                print "$name\t0\t-\n";
            }
        }
    }
}
close B;