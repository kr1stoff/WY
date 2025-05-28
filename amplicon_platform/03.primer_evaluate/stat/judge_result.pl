use strict;

my $merge = shift;
my $ref = shift;

open A,$merge || die $!;
while(<A>){
    chomp;
    next if /^species/;
    my @line = split /\t/,$_;
    my ($Left_Seq,$Right_Seq,$Amplicon_Name,$Chrom,$Amp_Start,$Amp_End,$Amplicon_Seq_add_before,$Amplicon_Seq,$Amplicon_Seq_add_after) = @line[29,30,31,32,33,34,47,48,49];
    my ($nei_f,$nei_r,$nei_amp_len,$f_len,$r_len,$nei_seq_add,$nei_seq) = @line[2,3,9,13,14,20,21];
    # print "$nei_seq\n";
    my $f_seq = substr($nei_seq,0,$f_len);
    my $r_seq = substr($nei_seq,-$r_len,$r_len);
    my $r_seq_rev = Complement_Reverse($r_seq);
    # print "$f_seq\t$nei_f\t$r_seq_rev\t$nei_r\n";

    my $f_seq = up($f_seq);
    my $nei_f = up($nei_f);
    my $r_seq_rev = up($r_seq_rev);
    my $nei_r = up($nei_r);
    my $nei_result;
    if ($f_seq eq $nei_f && $r_seq_rev eq $nei_r){
        $nei_result = "Nei_YES";
    }
    else{
        $nei_result = "Nei_NO";
    }

    my ($f,$block,$r);
    if ($Amplicon_Seq =~ /\[(\S+)\](\S+)\[(\S+)\]/){
        $f = $1;
        $block = $2;
        $r = $3;
    }
    $f = up($f);
    $block = up($block);
    $r = up($r);
    $Amplicon_Seq_add_before = up($Amplicon_Seq_add_before);
    $Amplicon_Seq_add_after = up($Amplicon_Seq_add_after);
    my $r_rev = Complement_Reverse($r);
    my $seq_merge = $Amplicon_Seq_add_before.$f.$block.$r.$Amplicon_Seq_add_after;
    my $start20 = $Amp_Start-20;
    my $end20 = $Amp_End+20;
    my $pos = "$Chrom:$start20-$end20";
    &getfa($pos);
    my $seq_faidx = `sed -n '2p' out.fa`;
    $seq_faidx =~ s/\n//;
    $seq_faidx = up($seq_faidx);

    `rm out.fa`;
    if ($f eq $Left_Seq && $r_rev eq $Right_Seq && $seq_merge eq $seq_faidx){
        print "$Amplicon_Name\t$nei_result\tYES\t";
    }
    else{
        print "$Amplicon_Name\t$nei_result\tNO\t";
    }
    print "$f\t$Left_Seq\t$r_rev\t$Right_Seq\t$seq_merge\t$seq_faidx\n";

}

sub getfa{
    #$_[0]：参考基因组
    #$_[1]：bed区间文件
    `samtools faidx $ref $_[0] | perl /home/lanlei/format_fasta.pl - > out.fa`;
}

sub Complement_Reverse{
        my $seq=shift;
        $seq=~tr/AGCTagct/TCGAtcga/;
        $seq=reverse($seq);
        return $seq;

}

sub up{
        my $seq=shift;
        $seq=~tr/agct/AGCT/;
        return $seq;

}