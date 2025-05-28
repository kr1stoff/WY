use strict;
use Digest::MD5 qw(md5_hex);

my $spe = shift;
my $inc = shift;
my $pri = shift;
my $taxid = shift;
my $fa = shift;
my $spe_cutoff = shift;
my $inc_cutoff = shift;

$spe_cutoff ||= 90;
$inc_cutoff ||= 85;

# 引物文件0base:0base

# 1	    ID	79a0b61de83883a1d598064d40317f53
# 2	    Left_Primer_Seq	TTTTTAAGGTTTAGAGCTAGAAATTTGC
# 3	    Right_Primer_Seq	ACCTGAACCTGCAAACAAATCTA
# 4	    Amplicon_Name	Streptococcus pyogenes-rsmD-5
# 5	    Chrom	NZ_LS483338.1
# 6	    Amplicon_Start	607240
# 7	    Amplicon_End	607439
# 8	    Amplicon_Length	200
# 9	    Amplified_Block_Start	607268
# 10	Amplified_Block_End	607416
# 11	Left_Primer_Tm	61
# 12	Right_Primer_Tm	62
# 13	Left_Primer_Length	28
# 14	Right_Primer_Length	23
# 15	Left_GC	29
# 16	Right_GC	39
# 17	Desired_Multiplex	True
# 18	Designed_Multiplex	1
# 19	Taxid	1314
# 20	Amplicon_Seq	TTTTAAGGTTTAGAGCTAGAAATTTGCTATAATAAATGTTATGAGAGTTGTATCAGGTGAATTTGGTGGGCGTCCTTTAAAGACGCTCGATGGAAAGATAACACGCCCTACTTCAGACAAAGTTAGAGGTGCGATTTTTAATATGATAGGACCTTACTTTAATGGTGGTCGTGTTTTAGATTTGTTTGCAGGTTCAGGT
# 21	1-specificity	100.000
# 22	2-specificity	91.815
# 23	inclusiveness	94.749
# 24	use	YES

open OUT,">out.primer.check.xls" or die $!;

print OUT "ID\tLeft_Seq\tRight_Seq\tAmplicon_Name\tChrom\tAmp_Start\tAmp_End\tAmp_Length\tAmplified_Block_Start\tAmplified_Block_End\tLeft_Tm\tRight_Tm\tLeft_Length\tRight_Length\tLeft_GC\tRight_GC\tDesired_Multiplex\tDesigned_Multiplex\tTaxid\tAmplicon_Seq\tPair_perecnt\tPair_total\tF_percent\tF_total\tR_percent\tR_total\tPrimer_perecnt\tPrimer_total\tGenome_num\tUse\n";

my (%spe,%inc);
open S,$spe || die $!;
while(<S>){
    chomp;
    my @line = split /\t/,$_;
    my $name = $line[0];
    $spe{$name}  = join "\t",@line[1..$#line];
}
close S;

open I,$inc || die $!;
while(<I>){
    chomp;
    my @line = split /\t/,$_;
    my $name = $line[0];
    $inc{$name}  = join "\t",@line[1..$#line];
}
close I;

open PRI,$pri || die $!;
while(<PRI>){
    chomp;
    next if $_ =~ /^#/;
    my @line = split/\t/,$_;
    my $pri_info;#引物设计数据
    my $length = (scalar @line) -1;

    if($length == 17){
        $pri_info = join("\t",@line[1..17]);
    }else{
        $pri_info = join("\t",@line[1..16],"-");
    }

    my $primer = join("-",$line[1],$line[2]);
    my $md5_primer = md5_hex($primer);
    my @name_tmps = split/-/,$line[3];
    my $name = join ('-',@name_tmps[1..$#name_tmps]);  ## 名字
    my $use; # 引物是否可用
    my $specificity = (split /\t/,$spe{$name})[0];
    my $inclusiveness = (split /\t/,$inc{$name})[0];
    # print "$specificity\t$inclusiveness\n";
    my ($chr,$amp_start,$amp_end) = @line[4,5,6];
    my $amp_start_0base = $amp_start - 1;
    my $amp_end_1base = $amp_end + 1;

    my $fa_id = "$line[4]:$amp_start_0base-$amp_end";
    my %fa;
    read_fasta($fa,\%fa);
    my $amp_seq = $fa{$fa_id};


    if ($specificity >= $spe_cutoff && $inclusiveness >= $inc_cutoff){
        $use = "YES";
    }
    else{
        $use = "NO";
    }

    print OUT "$md5_primer\t$pri_info\t$taxid\t$amp_seq\t$spe{$name}\t$inc{$name}\t$use\n";
}
close OUT;

sub read_fasta{
    my ($file,$bait_hp) = @_;
    open(IN,$file)||die("fail to open $file\n");
    $/=">";<IN>;$/="\n";
    while (<IN>) {
        my $title=$_;
        chomp $title;
        $/=">";
        my $seq=<IN>;
        chomp $seq;
        $seq =~ s/[\n\r]*//g;
        $/="\n";
        my @temp=split(/\s+/,$title);
        my $id = $temp[0];
        $$bait_hp{$id}=$seq;
    }
    close(IN);
}


