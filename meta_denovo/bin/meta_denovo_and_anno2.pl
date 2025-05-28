#!/usr/bin/perl -w 
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Getopt::Long;
use File::Basename qw(fileparse basename dirname);
use GACP qw(parse_config);
use Cwd qw(getcwd abs_path);

=head1 Description
    Genome assembly and anno pipeline.
    This pipeline is used to
    setp1. fastq file quality check ,Fastq filter and rm host
    step2. genome assembly and assembly result Assessment
    step3. Contig abundance and Contig taxonomy
    step4. Gene prediction and abundance
    step5. Gene function analysis
    setp6. Report analysis

=head1 Usage
    perl meta_denovo_and_anno.pl
    -fq_input <s>           The file list of raw fastq, format:name raw_fq1 raw_fq2; required if step1.
    -clean_input <s>        The file list of clean fastq, format:name clean_fq1 clean_fq2; required if step2.
    -workdir <s>            The work directory for analysis, default: "./".
    -step <s>               The step of the pipeline, defaul: 12345.
    -filter_software <s>    The  filter software (trimmomatic | fastp), default: fastp.
    -assembly_software <s>  The genome assembly software (spades | megahit), default: spades.
    -function_db <s>        The function database, (all | normal(default) | databese delimited by ",").
                            all: GO,NR2GO,KEGG,COG,eggNOG,CAZY,,VFDB,CARD,PHIbase,TCDB.
                            normal: GO,eggNOG2GO,KEGG,COG,eggNOG,CAZY,VFDB,CARD (default).
    -help                   The help message.

=head1 Exmple
    perl meta_denovo_and_anno.pl -fq_input fq.list -workdir ./ -step 123456 -filter_software fastp -assembly_software spades
    perl meta_denovo_and_anno.pl -clean_input clean.list -workdir ./ -step 123456
=cut


#################### Set global variables ####################
my ($fq_input,$clean_input,$genome_input,$bam_input,$pep_input,$workdir,$step,$Help);
# Config file

my ($filter_software,$assembly_software,$function_db);
my ($filter_dir,$assembly_dir,$contig_dir,$gene_dir,$function_dir,$report_dir);

GetOptions(
    "fq_input:s"=>\$fq_input,
    "clean_input:s"=>\$clean_input,
    "workdir:s"=>\$workdir,
    "step:s"=>\$step,
    "filter_software:s"=>\$filter_software,
    "assembly_software:s"=>\$assembly_software,
    "function_db:s"=>\$function_db,
    "help"=>\$Help
);

die `pod2text $0` if($Help || (! $fq_input && ! $clean_input));

$workdir ||= getcwd();
$workdir = abs_path($workdir);

$function_db ||= "GO,eggNOG2GO,KEGG,COG,eggNOG,CAZY,VFDB,CARD";
$step ||= "123456";
$filter_software ||= "fastp";
$assembly_software ||= "spades";
my $shell_dir="$workdir/00.shell";
`mkdir $shell_dir` unless (-d $shell_dir);

my $config_file = "$Bin/../config.txt";

# database
my $human_db = parse_config($config_file,"human_db");
my $nt = parse_config($config_file,"nt");
my $taxdump_db = parse_config($config_file,"taxdump_db");

# filter
my $fastqc = parse_config($config_file,"fastqc");
my $trimmomatic = parse_config($config_file,"trimmomatic");
my $fastp = parse_config($config_file,"fastp");

# assembly
my $spades = parse_config($config_file,"spades");
#my $megahit = parse_config($config_file,"megahit");

# contig
my $bowtie2 = parse_config($config_file,"bowtie2");
my $taxonkit = parse_config($config_file,"taxonkit");
# my $quast = parse_config($config_file,"quast"); 
my $checkM = parse_config($config_file,"checkM");
my $coverm = parse_config($config_file,"coverm");


# gene
#my $prokka = parse_config($config_file,"prokka");
my $bedtools = parse_config($config_file,"bedtools");

# func


my $time = `date`; chomp($time);

#################### End get parameters from command line ####################

my @shell_list;
if($step=~/1/){
    $filter_dir = "$workdir/01.Filter";
    `mkdir -p $filter_dir`;

    my $clean_list = "$workdir/clean.list";
    open IN1,"<$fq_input" || die $!;
    open LIST1,">$clean_list" || die $!;

    while(<IN1>){
        chomp;
        my ($name,$f1,$f2)=split /\s+/,$_;

        $name =~ s/^\s+|\s+$//g;
        my ($fq1,$fq2);

        my $dirname="$filter_dir/$name";
        `mkdir -p $dirname`;
        my $cmd = "cd $dirname \n";

        if ($f1 =~ /https/){
            $cmd .= "wget -O /sdbb/bioinfor/EARTH/RawData/$name.fq1.gz $f1 \n";
            $cmd .= "wget -O /sdbb/bioinfor/EARTH/RawData/$name.fq2.gz $f2 \n";
            $fq1 =  "/sdbb/bioinfor/EARTH/RawData/$name.fq1.gz";
            $fq2 =  "/sdbb/bioinfor/EARTH/RawData/$name.fq2.gz";
        }
        else{
            $fq1 = abs_path($f1);
            $fq2 = abs_path($f2);
        }

        # fastqc
        $cmd .= "ln -sf $fq1 $dirname/$name.raw.R1.fastq.gz\n";
        $cmd .= "ln -sf $fq2 $dirname/$name.raw.R2.fastq.gz\n";
        $cmd .= "$fastqc $dirname/$name.raw.R1.fastq.gz $dirname/$name.raw.R2.fastq.gz -t 4\n";
        if ($filter_software eq "trimmomatic"){
        # trimmomatic
            $cmd .= "java -jar $trimmomatic PE $fq1 $fq2 ";
            $cmd .= "$dirname/$name.clean.R1.fastq.gz $dirname/$name.unpaired.R1.fastq.gz $dirname/$name.clean.R2.fastq.gz $dirname/$name.unpaired.R2.fastq.gz ";
            $cmd .= "CROP:140 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 MINLEN:100 -threads 4\n";
        }
        if ($filter_software eq "fastp"){
        # fastp
            $cmd .="$fastp -i $fq1 -I $fq2 -q 15 -u 40 -n 5 -l 80 -w 4 -o $dirname/$name.clean.R1.fastq.gz -O $dirname/$name.clean.R2.fastq.gz 2> report.log\n";
        }
        # fastqc
        $cmd .= "$fastqc $dirname/$name.clean.R1.fastq.gz $dirname/$name.clean.R2.fastq.gz -t 4\n";

        # rmhost
        $cmd .= "$bowtie2 -p 10 -x $human_db -1 $dirname/$name.clean.R1.fastq.gz -2 $dirname/$name.clean.R2.fastq.gz -S $dirname/$name.sam --un-conc $dirname/$name.fq \n";

        # generate work shell
        generateShell("$shell_dir/S1_filter.$name.sh", $cmd, \@shell_list);
#        push @shell_list, "$shell_dir/S1_filter.$name.sh";
        print LIST1 "$name\t$dirname/$name.1.fq\t$dirname/$name.2.fq\n";
    }
    close IN1;
    close LIST1;
}

if($step=~/2/){
    $assembly_dir = "$workdir/02.Assembly";
    `mkdir -p $assembly_dir`;

    $clean_input = "$workdir/clean.list" if (! defined $clean_input);
    my $genome_list = "$workdir/genome.list";
    open IN2, "<$clean_input" || die $!;
    open LIST2,">$genome_list" || die $!;
    while(<IN2>){
        chomp;
        my ($name,$fq_clean1,$fq_clean2)=split /\s+/,$_;
        $fq_clean1 = abs_path($fq_clean1);
        $fq_clean2 = abs_path($fq_clean2);

        my $dirname="$assembly_dir/$name";
        `mkdir -p $dirname`;
        my $cmd = "cd $dirname \n";
        if ($assembly_software eq "spades"){
        # spades
            $cmd .= "$spades -t 20 --meta -1 $fq_clean1 -2 $fq_clean2 -o spades_output \n";
            $cmd .= "perl $Bin/assembly/spades_fa_format.pl $dirname/spades_output/scaffolds.fasta > $dirname/$name.genome.fa \n";
            $cmd .= "bowtie2-build --quiet $dirname/$name.genome.fa $dirname/$name.genome.fa";
        }
        # if ($assembly_software eq "megahit"){
        # # megahit
        #     $cmd .= "$SOAPdenovo2 all -s lib.list -K $soap_kmer -d 1 -D 1 -o kmer$soap_kmer -F >kmer$soap_kmer.log\n";
        #     $cmd .= "cp $dirname/kmer$soap_kmer.scafSeq $dirname/$name.genome.fa\n";
        # }

        generateShell("$shell_dir/S2_assembly.$name.sh", $cmd, \@shell_list);
        #push @shell_list, "$shell_dir/S2_assembly.$name.sh";
        print LIST2 "$name\t$dirname/$name.genome.fa\t$fq_clean1\t$fq_clean2\n";
    }
    close IN2;
    close LIST2;
}

if($step=~/3/){
    $contig_dir = "$workdir/03.Contig_abundance_and_taxonomy";
    `mkdir -p $contig_dir`;

    $genome_input = "$workdir/genome.list";
    open IN3, "<$genome_input" || die $!;

    $bam_input = "$workdir/bam.list";
    open LIST3,">$bam_input" || die $!;
    
    while(<IN3>){
        chomp;
        my ($name,$genome,$fq_clean1,$fq_clean2)=split /\s+/,$_;
        $genome = abs_path($genome);
        my $filename = basename($genome);
        
        # get the suffix of the genome fa for checkM.
        my @suffixlist=qw(fa fasta fna);
        my $suffix = (fileparse($genome,@suffixlist))[-1];

        my $dirname="$contig_dir/$name";
        `mkdir -p $dirname`;
        # quast
        my $abundance_dir = "$dirname/01.abundance";
        my $cmd_abundance = "mkdir -p $abundance_dir\n";
        $cmd_abundance .= "cd $abundance_dir\n";
        $cmd_abundance .= "$bowtie2 -p 10 -x $genome -1 $fq_clean1 -2 $fq_clean2 | samtools view --threads 10 -bS - > $abundance_dir/$name.bam\n";
        $cmd_abundance .= "samtools sort --threads 10 $abundance_dir/$name.bam -o $abundance_dir/$name.sort.bam\n";
        $cmd_abundance .= "samtools index $abundance_dir/$name.sort.bam\n";
        $cmd_abundance .= "samtools flagstat $abundance_dir/$name.sort.bam > $abundance_dir/$name.flagstat\n";
        $cmd_abundance .= "perl $Bin/gene/get_total_reads_stat.pl $abundance_dir/$name.flagstat > $abundance_dir/$name.total.reads \n";
        $cmd_abundance .= "$coverm filter --bam-files $abundance_dir/$name.sort.bam --output-bam-files $abundance_dir/$name.sort.filter.bam --threads 10 --min-read-percent-identity 90\n";
        $cmd_abundance .= "samtools index $abundance_dir/$name.sort.filter.bam \n";
        $cmd_abundance .= "$coverm contig --bam-files $abundance_dir/$name.sort.filter.bam -m trimmed_mean --threads 10 --output-file $abundance_dir/$name.sort.filter.tpmean\n";
        $cmd_abundance .= "$bedtools genomecov -bga -pc -ibam $abundance_dir/$name.sort.filter.bam > $abundance_dir/$name.sort.filter.cov\n";
        $cmd_abundance .= "perl $Bin/contig/filter_contig_cov.pl $abundance_dir/$name.sort.filter.cov > $abundance_dir/$name.sort.filter.cov.contig\n";
        $cmd_abundance .= "perl $Bin/contig/fishInWinter.pl $abundance_dir/$name.sort.filter.cov.contig $abundance_dir/$name.sort.filter.tpmean > $abundance_dir/$name.sort.filter.cov.contig.abundance\n";
        
        generateShell("$shell_dir/S3_01_abundance.$name.sh", $cmd_abundance, \@shell_list);
        #push @shell_list, "$shell_dir/S3_01_abundance.$name.sh";


        ###################################################################################################################
        # taxonomy
        my $taxonomy_dir = "$dirname/02.taxonomy";
        my $cmd_taxonomy = "mkdir -p $taxonomy_dir\n";

        $cmd_taxonomy .= "cd $taxonomy_dir\n";
        $cmd_taxonomy .= "export BLASTDB=/sdbb/bioinfor/Database/nt/2021-07-21/ \n";
        $cmd_taxonomy .= "blastn -db $nt/nt -query $genome -evalue 1e-10 -num_threads 20 ";
        $cmd_taxonomy .= "-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus staxid sscinames scomnames' -out $taxonomy_dir/$name.blastn.out\n";
        $cmd_taxonomy .= "awk '\$15>50' $taxonomy_dir/$name.blastn.out | perl $Bin/contig/get_blast_top_n.pl - 5 > $taxonomy_dir/$name.blastn.out.top5 \n";
        $cmd_taxonomy .= "cut -f 1,18 $taxonomy_dir/$name.blastn.out.top5 | perl $Bin/contig/get_lca_input.pl > $taxonomy_dir/$name.contig2taxids.txt \n";
        $cmd_taxonomy .= "cut -f 2 $taxonomy_dir/$name.contig2taxids.txt | $taxonkit lca --data-dir $taxdump_db  > $taxonomy_dir/contig.tmp \n";
        $cmd_taxonomy .= "paste $taxonomy_dir/$name.contig2taxids.txt $taxonomy_dir/contig.tmp |cut -f 1,4 > $taxonomy_dir/$name.contig.taxids.lca \n";
        $cmd_taxonomy .= "rm $taxonomy_dir/contig.tmp \n";
        generateShell("$shell_dir/S3_02_taxonomy.$name.sh", $cmd_taxonomy, \@shell_list);
        #push @shell_list, "$shell_dir/S3_02_taxonomy.$name.sh";

        print LIST3 ("$name\t$genome\t$abundance_dir/$name.sort.filter.bam\t$abundance_dir/$name.total.reads\t$fq_clean1\t$fq_clean2\n")

    }
    close IN3;
}

if($step=~/4/){
    $gene_dir = "$workdir/04.Gene";
    `mkdir -p $gene_dir`;

    $genome_input = "$workdir/bam.list";
    open IN4, "<$genome_input" || die $!;
    
    my $pep_list = "$workdir/pep.list";
    open LIST4,">$pep_list" || die $!;

    while(<IN4>){
        chomp;
        my ($name,$genome,$bam,$total_reads,$fq_clean1,$fq_clean2)=split /\s+/,$_;
        $genome = abs_path($genome);
        my $filename = basename($genome);

        my $dirname="$gene_dir/$name";
        `mkdir -p $dirname`;
        # gene
        my $gene_dir = "$dirname/01.gene_prediction";
        my $cmd_gene = "mkdir -p $gene_dir\n";
        $cmd_gene .= "cd $gene_dir\n";
        $cmd_gene .= "source /home/lanlei/miniconda3/bin/activate /home/lanlei/miniconda3/envs/prokka \n";
        #$cmd_gene .= "source /sdbb/bioinfor/zhuzi/soft/anaconda3/bin/activate  single_gene_assemble \n";
        $cmd_gene .= "prokka $genome --metagenome --prefix $name --outdir $gene_dir/$name --force --quiet\n";
        $cmd_gene .= "source /home/lanlei/miniconda3/bin/deactivate \n";
        $cmd_gene .= "cp $gene_dir/$name/$name.faa $dirname/$name.faa\n";
        $cmd_gene .= "cp $gene_dir/$name/$name.gff $dirname/$name.gff\n";

        generateShell("$shell_dir/S4_01_gene.$name.sh", $cmd_gene, \@shell_list);
        #push @shell_list, "$shell_dir/S4_01_gene.$name.sh";

        # gene_abun
        my $gene_abun_dir = "$dirname/02.gene_abun";
        my $cmd_gene_abun = "mkdir -p $gene_abun_dir\n";
        $cmd_gene_abun .= "cd $gene_abun_dir\n";
        $cmd_gene_abun .= "cat $dirname/$name.gff | grep CDS | cut -f 1,4,5,9 | awk -F ';' '{print \$1}' | sed 's/ID=//' > $gene_abun_dir/$name.gene.bed \n";
        $cmd_gene_abun .= "bedtools multicov -bams $bam -bed $gene_abun_dir/$name.gene.bed > $gene_abun_dir/$name.gene.count \n";
        $cmd_gene_abun .= "perl $Bin/gene/cal_RPKM.pl $total_reads $gene_abun_dir/$name.gene.count |  cut -f 4,8 > $gene_abun_dir/$name.rpkm \n";
        generateShell("$shell_dir/S4_02_gene_abun.$name.sh", $cmd_gene_abun, \@shell_list);
        #push @shell_list, "$shell_dir/S4_02_gene_abun.$name.sh";

        print LIST4 "$name\t$dirname/$name.faa\n";
    }
    close IN4;
    close LIST4;
}

if($step=~/5/){
    $function_dir = "$workdir/05.function";
    `mkdir -p $function_dir`;

    $pep_input = "$workdir/pep.list";
    open IN5, "<$pep_input" || die $!;

    while(<IN5>){
        chomp;
        my ($name,$pep)=split /\s+/,$_;
        $pep = abs_path($pep);

        my $dirname="$function_dir/$name";
        if ($function_db eq "normal"){
            $function_db = "GO,eggNOG2GO,KEGG,COG,eggNOG,CAZY,VFDB,CARD";
        }
        elsif ($function_db eq "all"){
            $function_db = "GO,NR2GO,KEGG,COG,eggNOG,CAZY,VFDB,CARD,PHIbase,TCDB";
        }
        my @func = split /,/,$function_db;
        my $all_func_db = "-" . join (" -",@func);

        my $cmd_func = "mkdir -p $dirname \n";
        $cmd_func .= "cd $dirname \n";
        $cmd_func .= "perl $Bin/function_anno_auto.pl $all_func_db $pep\n";
        $cmd_func .= "bash $function_dir/$name/STEP01_fun_ann_work.sh\n";
        $cmd_func .= "bash $function_dir/$name/STEP02_fun_ann_stat.sh\n";
        #$cmd_func .= "sh $function_dir/$name/STEP03_delete_tmp_files.sh\n";
        generateShell("$shell_dir/S5_function.$name.sh", $cmd_func ,\@shell_list);
    }
    close IN5;
}

if($step=~/6/){
    $report_dir = "$workdir/report";
    `mkdir -p $report_dir`;

    $pep_input = "$workdir/pep.list";
    open IN6, "<$pep_input" || die $!;

    while(<IN6>){
        chomp;
        my ($name,$pep)=split /\s+/,$_;
        $pep = abs_path($pep);

        my $dirname = "$report_dir/$name";
        my $cmd_report;
        
        $cmd_report .= "mkdir -p $dirname/01.Filter/ $dirname/02.Assembly/ $dirname/05.function/ $dirname/report/\n";
        $cmd_report .= "mkdir -p $dirname/03.Contig_abundance_and_taxonomy/01.abundance $dirname/03.Contig_abundance_and_taxonomy/02.taxonomy \n";
        $cmd_report .= "mkdir -p $dirname/04.Gene/01.gene_prediction $dirname/04.Gene/02.gene_abun \n";
        $cmd_report .= "cd $dirname/05.function && mkdir CARD COG eggNOG GO KEGG \n";

        $cmd_report .= "cp $workdir/01.Filter/$name/*.html $dirname/01.Filter/\n";
        $cmd_report .= "cp $workdir/03.Contig_abundance_and_taxonomy/$name/01.abundance/*.sort.filter.cov.contig.abundance $dirname/03.Contig_abundance_and_taxonomy/01.abundance/${name}_contig_aubn.xls \n";
        $cmd_report .= "cp $workdir/03.Contig_abundance_and_taxonomy/$name/02.taxonomy/*.contig.taxids.lca $dirname/03.Contig_abundance_and_taxonomy/02.taxonomy/${name}_Contig_taxon.xls \n";
        $cmd_report .= "cp $workdir/04.Gene/$name/02.gene_abun/*rpkm  $dirname/04.Gene/02.gene_abun/${name}_gene_abun.xls \n";
        $cmd_report .= "cp $workdir/05.function/$name/CARD/output.txt  $dirname/05.function/CARD/${name}_card.xls \n";
        $cmd_report .= "cp $workdir/05.function/$name/COG/*.faa.cog.xls $workdir/05.function/$name/COG/*result.COG.png  $dirname/05.function/COG/ \n";
        $cmd_report .= "cp $workdir/05.function/$name/eggNOG/out.emapper.annotations $dirname/05.function/eggNOG/ \n";
        # go
        $cmd_report .= "cp $workdir/05.function/$name/GO_eggNOG/*.faa.GO2Gene.xls $dirname/05.function/GO/${name}_GO2Gene.xls \n";
        $cmd_report .= "cp $workdir/05.function/$name/GO_eggNOG/*.faa.GO.png $dirname/05.function/GO/${name}_GO.png \n";
        $cmd_report .= "cp $workdir/05.function/$name/KEGG/*.faa.kegg.xls  $dirname/05.function/KEGG/${name}_kegg.xls \n";
        $cmd_report .= "cp /sdbb/bioinfor/lanlei/meta_denovo/9.report/meta_denovo/meta_denovo_report.pl  $dirname/report/ \n";
        $cmd_report .= "cp -R /sdbb/bioinfor/lanlei/meta_denovo/9.report/meta_denovo/src $dirname/report/ \n";
        $cmd_report .= "cd  $dirname/report/ \n";
        $cmd_report .= "perl meta_denovo_report.pl $name \n";
        
        generateShell("$shell_dir/S6_report.$name.sh", $cmd_report ,\@shell_list);

    }
    close IN6;
}

open MAIN,">$shell_dir/main.sh" || die $!;
print MAIN join (" && \\\n",@shell_list) . "\n";
close MAIN;


sub generateShell{
    my ($output_shell, $content, $outshell ,$finish_string) = @_;
    unlink glob "$output_shell.*";
    $finish_string ||= "Still_waters_run_deep";
    chomp $content;
    open OUT,">$output_shell" or die "Cannot open file $output_shell:$!";
    print OUT "#!/bin/bash\n";
    print OUT "echo ==========start at : `date` ==========\n";
    print OUT "set -e \n";
    print OUT "$content && \\\n";
    print OUT "echo ==========end at : `date` ========== && \\\n";
    print OUT "echo $finish_string 1>&2 && \\\n";
    print OUT "echo $finish_string > $output_shell.sign\n";
    close OUT;
    my $run = "bash $output_shell 1>$output_shell.e 2>$output_shell.o";
    push @{$outshell}, $run;
}
