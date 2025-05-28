use strict;

my $csv = shift;

my @all_file;
open A,$csv || die $!;
while(<A>){
    chomp;
    my @a = split /\t/,$_;
    push @all_file,$a[1];
    if (exists $a[2]){
        push @all_file,$a[2];
    }
}
close A;

foreach my $file (@all_file){
    if(-e $file){
        if (-z $file){
            print "$file is null. \n";
        }
        else{
            if (-d $file){
                print "$file is a directory, not a file. \n";
            }
        }
    }
    else{
	    print "$file do not exist. \n";
    } 
}
