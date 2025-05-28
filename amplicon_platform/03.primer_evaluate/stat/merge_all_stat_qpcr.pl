use strict;
use Digest::MD5 qw(md5_hex);

my $spe = shift;
my $inc = shift;
my $exp = shift;
my $taxid = shift;
my $fa = shift;
my $spe_cutoff = shift;
my $inc_cutoff = shift;

$spe_cutoff ||= 90;
$inc_cutoff ||= 85;

# 1	#	1
# 2	Fwd Start	459
# 3	Fwd Stop	476
# 4	Fwd Length	18
# 5	Fwd Tm	59
# 6	Fwd %GC	56
# 7	Fwd Seq	TCGGAGGCGAGAAATCCA
# 8	Rev Start	574
# 9	Rev Stop	552
# 10	Rev Length	23
# 11	Rev Tm	60
# 12	Rev %GC	43
# 13	Rev Seq	AAAAGCTTACCGGCTCATCAATC
# 14	Probe Start	487
# 15	Probe Stop	506
# 16	Probe Length	20
# 17	Probe Tm	69
# 18	Probe %GC	65
# 19	Probe Seq	CCATGCCGCCGCCATTACGT
# 20	Amp Tm	85
# 21	Amp %GC	62
# 22	Amp Ta	62
# 23	Amp Len	116
# 24	Penalty	335
# 25	Chromosome	NZ_CP008782.1
# 26	Start	2040086
# 27	End	2040692
# 28	Name	Burkholderia pseudomallei-group_16135-1

open OUT,">qPCR.primer.check.xls" or die $!;
open TNGS,">2tNGS.bed" or die $!;

print OUT "ID\tLeft_Seq\tRight_Seq\tProbe_Seq\tAmplicon_Name\tChrom\tAmp_Start\tAmp_End\tAmp_Length\tLeft_Tm\tRight_Tm\tProbe_Tm\tLeft_Length\tRight_Length\tProbe_Length\tLeft_GC\tRight_GC\tProbe_GC\tTaxid\tAmplicon_Seq\tPair_perecnt\tPair_total\tF_percent\tF_total\tR_percent\tR_total\tPer_same\tProbe_Total\tPrimer_perecnt\tPrimer_total\tProbe_percent\tProbe_total\tGenome_num\tUse\n";

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

open PRI,$exp || die $!;
while(<PRI>){
    chomp;
    next if $_ =~ /^#/;
    next if $_ =~ /^Index/;
    my @line = split/\t/,$_;
    my $amp_start = $line[25] + $line[1]; # 此处不一样
    my $amp_start_0base = $amp_start - 1;
    my $amp_end = $line[25] + $line[7]; # 此处不一样
    my $exp_info = join("\t",$line[6],$line[12],$line[18],$line[27],$line[24],$amp_start,$amp_end,$line[22],$line[4],$line[10],$line[16],$line[3],$line[9],$line[15],$line[5],$line[11],$line[17]);
    my $primer_probe = join("-",$line[6],$line[12],$line[18]);
    my $md5_primer_probe = md5_hex($primer_probe);
    my @name_tmps = split/-/,$line[27];
    my $name = join ('-',@name_tmps[1..$#name_tmps]);  ## 名字
    my $use; # 引物是否可用
    my $specificity = (split /\t/,$spe{$name})[0];
    my ($inclusiveness_primer,$inclusiveness_probe) = (split /\t/,$inc{$name})[0,2];
    # print "$specificity\t$inclusiveness_primer\t$inclusiveness_probe\n";
    my $fa_id = "$line[24]:$amp_start_0base-$amp_end";
    my %fa;
    read_fasta($fa,\%fa);
    my $amp_seq = $fa{$fa_id};


    if ($specificity >= $spe_cutoff && $inclusiveness_primer >= $inc_cutoff && $inclusiveness_probe >= $inc_cutoff){
        $use = "YES";
        print TNGS "$line[24]\t$amp_start\t$amp_end\t$line[27]\t1\tB\tW\n";
    }
    else{
        $use = "NO";
    }

    print OUT "$md5_primer_probe\t$exp_info\t$taxid\t$amp_seq\t$spe{$name}\t$inc{$name}\t$use\n";
}

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





