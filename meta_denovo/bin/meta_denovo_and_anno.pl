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
    -help                   The help message.

=head1 Exmple
    perl meta_denovo_and_anno.pl -fq_input fq.list -workdir ./ -step 123456 -filter_software fastp -assembly_software spades
    perl meta_denovo_and_anno.pl -clean_input clean.list -workdir ./ -step 123456
=cut


#################### Set global variables ####################
my ($fq_input,$clean_input,$genome_input,$bam_input,$pep_input,$workdir,$step,$Help);
# Config file

my ($filter_software,$assembly_software,$task_id);
my ($filter_dir,$assembly_dir,$contig_dir,$gene_dir,$function_dir,$Upload_dir,$log_dir);

GetOptions(
    "fq_input:s"=>\$fq_input,
    "clean_input:s"=>\$clean_input,
    "workdir:s"=>\$workdir,
    "step:s"=>\$step,
    "filter_software:s"=>\$filter_software,
    "assembly_software:s"=>\$assembly_software,
    "task_id:s"=>\$task_id,
    "help"=>\$Help
);


die `pod2text $0` if($Help || (! $fq_input && ! $clean_input));

$workdir ||= getcwd();
$workdir = abs_path($workdir);

$step ||= "123456";
$filter_software ||= "fastp";
$assembly_software ||= "spades";
my $shell_dir="$workdir/00.shell";
`mkdir $shell_dir` unless (-d $shell_dir);

$log_dir = "$workdir/log";
`mkdir -p $log_dir`;

my $config_file = "$Bin/../config.txt";

# base
my $python = parse_config($config_file,"python");
my $perl = parse_config($config_file,"perl");
my $Rscript = parse_config($config_file,"Rscript");

# database
my $human_db = parse_config($config_file,"human_db");
my $blast_db = parse_config($config_file,"blast_db");
my $nt = parse_config($config_file,"nt");
my $taxdump_db = parse_config($config_file,"taxdump_db");

# filter
my $fastqc = parse_config($config_file,"fastqc");
# my $trimmomatic = parse_config($config_file,"trimmomatic");
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
my $bedtools = parse_config($config_file,"bedtools");

# func
my $vfdb_dir = parse_config($config_file,"vfdb_dir");

# file
my $activate = parse_config($config_file,"activate");
my $AnalysisResults = parse_config($config_file,"AnalysisResults");


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
        my ($name,$fq1,$fq2)=split /\s+/,$_;

        $name =~ s/^\s+|\s+$//g;

        my $filter_each_dir="$filter_dir/$name";
        `mkdir -p $filter_each_dir`;
        my $cmd = "cd $filter_each_dir \n";

        $fq1 = abs_path($fq1);
        $fq2 = abs_path($fq2);


        # fastqc
        my ($fastqc_raw1,$fastqc_raw2);

        if ($fq1 =~ /gz/){
            $cmd .= "ln -sf $fq1 $filter_each_dir/$name.raw.R1.fastq.gz\n";
            $fastqc_raw1 = "$filter_each_dir/$name.raw.R1.fastq.gz";
            $cmd .= "$fastqc $filter_each_dir/$name.raw.R1.fastq.gz -o $filter_each_dir --extract --quiet -t 20\n";

            $cmd .= "ln -sf $fq2 $filter_each_dir/$name.raw.R2.fastq.gz\n";
            $fastqc_raw2 = "$filter_each_dir/$name.raw.R2.fastq.gz";
            $cmd .= "$fastqc $filter_each_dir/$name.raw.R2.fastq.gz -o $filter_each_dir --extract --quiet -t 20\n";
        }
        else{
            $cmd .= "ln -sf $fq1 $filter_each_dir/$name.raw.R1.fastq\n";
            $fastqc_raw1 = "$filter_each_dir/$name.raw.R1.fastq";
            $cmd .= "$fastqc $filter_each_dir/$name.raw.R1.fastq -o $filter_each_dir --extract --quiet -t 20\n";

            $cmd .= "ln -sf $fq2 $filter_each_dir/$name.raw.R2.fastq\n";
            $fastqc_raw2 = "$filter_each_dir/$name.raw.R2.fastq";
            $cmd .= "$fastqc $filter_each_dir/$name.raw.R2.fastq -o $filter_each_dir --extract --quiet -t 20\n";
        }

        if ($filter_software eq "fastp"){
        # fastp
            my $adapter = 'AGATCGGAAGAGC';
            $cmd .="$fastp -a $adapter -i $fq1 -I $fq2 -q 15 -u 40 -n 5 -l 30 -w 20 -j $filter_each_dir/$name.clean.fq.stat.json -o $filter_each_dir/$name.clean.R1.fastq.gz -O $filter_each_dir/$name.clean.R2.fastq.gz 2> fq_qc.log\n";
            $cmd .="$python $Bin/stat/parse_fastp.json.py -id $name -i $filter_each_dir/$name.clean.fq.stat.json -m PE \n";
        }   
        # fastqc
        $cmd .= "$fastqc $filter_each_dir/$name.clean.R1.fastq.gz $filter_each_dir/$name.clean.R2.fastq.gz -o $filter_each_dir --extract --quiet -t 20\n";

        # rmhost
        $cmd .= "$bowtie2 -p 20 -x $human_db -1 $filter_each_dir/$name.clean.R1.fastq.gz -2 $filter_each_dir/$name.clean.R2.fastq.gz -S $filter_each_dir/$name.sam --un-conc $filter_each_dir/$name.fq \n";

        # generate work shell
        generateShell("$shell_dir/S1_filter.$name.sh", $cmd, \@shell_list);
#        push @shell_list, "$shell_dir/S1_filter.$name.sh";
        print LIST1 "$name\t$filter_each_dir/$name.1.fq\t$filter_each_dir/$name.2.fq\n";
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

        my $assembly_each_dir="$assembly_dir/$name";
        `mkdir -p $assembly_each_dir`;
        my $cmd = "cd $assembly_each_dir \n";
        if ($assembly_software eq "spades"){
        # spades
            $cmd .= "$spades -t 40 --meta -1 $fq_clean1 -2 $fq_clean2 -o spades_output \n";
            $cmd .= "$perl $Bin/assembly/spades_fa_format.pl $assembly_each_dir/spades_output/scaffolds.fasta > $assembly_each_dir/$name.genome.fa \n";
            $cmd .= "$python $Bin/stat/ass_fa_stat.py -fa $assembly_each_dir/$name.genome.fa -p ${name} \n";
            $cmd .= "bowtie2-build --quiet $assembly_each_dir/$name.genome.fa $assembly_each_dir/$name.genome.fa \n";
        }
        # if ($assembly_software eq "megahit"){
        # # megahit
        #     $cmd .= "$SOAPdenovo2 all -s lib.list -K $soap_kmer -d 1 -D 1 -o kmer$soap_kmer -F >kmer$soap_kmer.log\n";
        #     $cmd .= "cp $assembly_each_dir/kmer$soap_kmer.scafSeq $assembly_each_dir/$name.genome.fa\n";
        # }

        # $cmd .= "$perl $Bin/stat/seq_n50.pl $assembly_each_dir/$name.genome.fa > $assembly_each_dir/$name.contig.stat \n";
        $cmd .= "$perl $Bin/stat/seq_len_freq.v4.pl -o $assembly_each_dir $assembly_each_dir/$name.genome.fa \n";

        generateShell("$shell_dir/S2_assembly.$name.sh", $cmd, \@shell_list);
        #push @shell_list, "$shell_dir/S2_assembly.$name.sh";
        print LIST2 "$name\t$assembly_each_dir/$name.genome.fa\t$fq_clean1\t$fq_clean2\n";
    }
    close IN2;
    close LIST2;
}

if($step=~/3/){
    $contig_dir = "$workdir/03.Contig";
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

        my $contig_each_dir="$contig_dir/$name";
        `mkdir -p $contig_each_dir`;
        # quast
        my $abundance_dir = "$contig_each_dir/01.abundance";
        my $cmd_abundance = "mkdir -p $abundance_dir\n";
        $cmd_abundance .= "cd $abundance_dir\n";
        $cmd_abundance .= "$bowtie2 -p 40 -x $genome -1 $fq_clean1 -2 $fq_clean2 | samtools view --threads 10 -bS - > $abundance_dir/$name.bam\n";
        $cmd_abundance .= "samtools sort --threads 40 $abundance_dir/$name.bam -o $abundance_dir/$name.sort.bam\n";
        $cmd_abundance .= "samtools index $abundance_dir/$name.sort.bam\n";
        $cmd_abundance .= "samtools flagstat $abundance_dir/$name.sort.bam > $abundance_dir/$name.flagstat\n";
        $cmd_abundance .= "$perl $Bin/gene/get_total_reads_stat.pl $abundance_dir/$name.flagstat > $abundance_dir/$name.total.reads \n";
        $cmd_abundance .= "$coverm filter --bam-files $abundance_dir/$name.sort.bam --output-bam-files $abundance_dir/$name.sort.filter.bam --threads 10 --min-read-percent-identity 90\n";
        $cmd_abundance .= "samtools index $abundance_dir/$name.sort.filter.bam \n";
        $cmd_abundance .= "$coverm contig --bam-files $abundance_dir/$name.sort.filter.bam -m trimmed_mean --threads 10 --output-file $abundance_dir/$name.sort.filter.tpmean\n";
        $cmd_abundance .= "$bedtools genomecov -bga -pc -ibam $abundance_dir/$name.sort.filter.bam > $abundance_dir/$name.sort.filter.cov\n";
        $cmd_abundance .= "$perl $Bin/contig/filter_contig_cov.pl $abundance_dir/$name.sort.filter.cov > $abundance_dir/$name.sort.filter.cov.contig\n";
        $cmd_abundance .= "$perl $Bin/contig/fishInWinter.pl $abundance_dir/$name.sort.filter.cov.contig $abundance_dir/$name.sort.filter.tpmean > $abundance_dir/$name.sort.filter.cov.contig.abundance\n";
        
        generateShell("$shell_dir/S3_01_abundance.$name.sh", $cmd_abundance, \@shell_list);
        #push @shell_list, "$shell_dir/S3_01_abundance.$name.sh";


        ###################################################################################################################
        # taxonomy
        my $taxonomy_dir = "$contig_each_dir/02.taxonomy";
        my $cmd_taxonomy = "mkdir -p $taxonomy_dir\n";

        $cmd_taxonomy .= "cd $taxonomy_dir\n";
        $cmd_taxonomy .= "export BLASTDB=$blast_db \n";
        $cmd_taxonomy .= "blastn -db $nt/nt -query $genome -evalue 1e-10 -num_threads 40 ";
        $cmd_taxonomy .= "-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus staxid sscinames scomnames' -out $taxonomy_dir/$name.blastn.out\n";
        $cmd_taxonomy .= "awk '\$15>30' $taxonomy_dir/$name.blastn.out | $perl $Bin/contig/get_blast_top_n.pl - 5 > $taxonomy_dir/$name.blastn.out.top5 \n";
        $cmd_taxonomy .= "cut -f 1,18 $taxonomy_dir/$name.blastn.out.top5 | $perl $Bin/contig/get_lca_input.pl > $taxonomy_dir/$name.contig2taxids.txt \n";
        $cmd_taxonomy .= "cut -f 2 $taxonomy_dir/$name.contig2taxids.txt | $taxonkit lca --data-dir $taxdump_db  > $taxonomy_dir/contig.tmp \n";
        $cmd_taxonomy .= "paste $taxonomy_dir/$name.contig2taxids.txt $taxonomy_dir/contig.tmp |cut -f 1,4 > $taxonomy_dir/$name.contig.taxids.lca \n";
        $cmd_taxonomy .= "cut -f 2 $taxonomy_dir/$name.contig.taxids.lca | taxonkit lineage --data-dir $taxdump_db - | taxonkit reformat -f \"{k};{p};{c};{o};{f};{g};{s}\" -F -P --data-dir $taxdump_db > $taxonomy_dir/$name.contig.taxids.lineage.tmp \n";
        $cmd_taxonomy .= "paste $taxonomy_dir/$name.contig.taxids.lca $taxonomy_dir/$name.contig.taxids.lineage.tmp | cut -f 1,5 > $taxonomy_dir/$name.contig.taxids.lineage.txt \n";
        $cmd_taxonomy .= "$perl $Bin/contig/combine_abundance_and_lineage.pl $taxonomy_dir/$name.contig.taxids.lineage.txt $abundance_dir/$name.sort.filter.cov.contig.abundance >  $contig_each_dir/${name}_Contig_abun_and_lineage.txt \n";
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
        my $gene_each_dir="$gene_dir/$name";
        `mkdir -p $gene_each_dir`;
        # gene
        my $gene_pre_dir = "$gene_each_dir/01.gene_prediction";
        my $cmd_gene = "mkdir -p $gene_pre_dir\n";
        $cmd_gene .= "cd $gene_pre_dir\n";
        $cmd_gene .= "source $activate denovo \n";
        $cmd_gene .= "prokka $genome --metagenome --prefix $name --outdir $gene_pre_dir/$name --force --quiet \n";
        $cmd_gene .= "source $activate base_env\n";
        $cmd_gene .= "cp $gene_pre_dir/$name/$name.ffn $gene_each_dir/$name.ffn \n";
        $cmd_gene .= "cp $gene_pre_dir/$name/$name.faa $gene_each_dir/$name.faa \n";
        $cmd_gene .= "cp $gene_pre_dir/$name/$name.gff $gene_each_dir/$name.gff \n";

        $cmd_gene .= "grep -v \"#\" $gene_each_dir/$name.gff | grep 'CDS' | cut -f 1,4,5,9 | awk -F ';' '{print \$1}'|sed 's/ID=//g'  > $gene_pre_dir/$name.gene.bed \n";
        $cmd_gene .= "bedtools coverage -a $gene_pre_dir/$name.gene.bed -b $bam -mean > $gene_pre_dir/$name.gene.bed.depth \n";
        $cmd_gene .= "bedtools coverage -a $gene_pre_dir/$name.gene.bed -b $bam > $gene_pre_dir/$name.gene.bed.cov \n";
        $cmd_gene .= "paste $gene_pre_dir/$name.gene.bed.depth  $gene_pre_dir/$name.gene.bed.cov | cut -f 4,5,13 > $gene_pre_dir/$name.gene.bed.mean_and_cov \n";

        $cmd_gene .= "$perl $Bin/stat/seq_len_freq.v4.pl -o $gene_pre_dir $gene_each_dir/$name.ffn \n";

        generateShell("$shell_dir/S4_01_gene.$name.sh", $cmd_gene, \@shell_list);
        #push @shell_list, "$shell_dir/S4_01_gene.$name.sh";

        # gene_abun
        my $gene_abun_dir = "$gene_each_dir/02.gene_abun";
        my $cmd_gene_abun = "mkdir -p $gene_abun_dir\n";
        $cmd_gene_abun .= "cd $gene_abun_dir\n";
        # $cmd_gene_abun .= "cat $gene_each_dir/$name.gff | grep 'CDS' | cut -f 1,4,5,9 | awk -F ';' '{print \$1}' | sed 's/ID=//' > $gene_abun_dir/$name.gene.bed \n";
        $cmd_gene_abun .= "bedtools multicov -bams $bam -bed $gene_pre_dir/$name.gene.bed > $gene_abun_dir/$name.gene.count \n";
        $cmd_gene_abun .= "$perl $Bin/gene/cal_RPKM.pl $total_reads $gene_abun_dir/$name.gene.count |  cut -f 4,6,8 > $gene_abun_dir/$name.rpkm \n";
        generateShell("$shell_dir/S4_02_gene_abun.$name.sh", $cmd_gene_abun, \@shell_list);
        #push @shell_list, "$shell_dir/S4_02_gene_abun.$name.sh";

        print LIST4 "$name\t$gene_each_dir/$name.faa\n";
    }
    close IN4;
    close LIST4;
}

if($step=~/5/){
    
    $function_dir = "$workdir/05.Function";
    `mkdir -p $function_dir`;

    $pep_input = "$workdir/pep.list";
    open IN5, "<$pep_input" || die $!;

    while(<IN5>){
        chomp;
        my ($name,$pep)=split /\s+/,$_;
        $pep = abs_path($pep);

        my $function_each_dir="$function_dir/$name";
        # card
        my $cmd_func = "mkdir -p $function_each_dir/CARD \n";
        $cmd_func .= "source $activate rgi \n";
        $cmd_func .= "rgi main --input_type protein --clean --input_sequence $pep --output_file $function_each_dir/CARD/rgi_result \n";
        $cmd_func .= "source $activate base_env \n";
        
        #20220701
        $cmd_func .="cut -f 1 $function_each_dir/CARD/rgi_result.txt |cut -f1 -d ' '|sed 's/ORF_ID/Contig/g' > $function_each_dir/CARD/contig.id \n";
        $cmd_func .="cut -f 8,10,11,15,16,21  $function_each_dir/CARD/rgi_result.txt > $function_each_dir/CARD/cut_drug.txt\n";
        $cmd_func .="paste  $function_each_dir/CARD/contig.id  $function_each_dir/CARD/cut_drug.txt > $function_each_dir/CARD/cut_id_drug.txt\n";
        $cmd_func .="python3 $Bin/stat/parse.rgi.py $function_each_dir/CARD/cut_id_drug.txt $gene_dir/$name/01.gene_prediction/$name.gene.bed.mean_and_cov \n";
        # $cmd_func .= "python3 $Bin/stat/rgi_json_explain.py $function_each_dir/CARD/rgi_result.json $function_each_dir/CARD/card_anno.out $gene_dir/$name/01.gene_prediction/$name.gene.bed.mean_and_cov \n";
        $cmd_func .= "mv  $function_each_dir/CARD/detail_drug_resistance.txt $function_each_dir/CARD/card_anno.txt\n";
        $cmd_func .= "mv  $function_each_dir/CARD/detail_drug_resistance.xlsx $function_each_dir/CARD/card_anno.xlsx\n";

        # VFDB
        $cmd_func .= "mkdir -p $function_each_dir/VFDB \n";
        $cmd_func .= "blastp -num_threads 60 -evalue 1e-10 -query $pep -db $vfdb_dir/VFDB -outfmt \"6 qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp\" -out $function_each_dir/VFDB/VFDB_anno.result \n";
        $cmd_func .= "$python $Bin/function/VFDB_filter_blastp.py $function_each_dir/VFDB/VFDB_anno.result $vfdb_dir/clean_VFDB_setB_des.txt $gene_dir/$name/01.gene_prediction/$name.gene.bed.mean_and_cov $function_each_dir/VFDB \n";
        $cmd_func .= "$python $Bin/function/get_balst_top.py -i $function_each_dir/VFDB/show_virulence_gene.txt -n 1 -o $function_each_dir/VFDB/top1_VFDB_anno.txt \n";

        $cmd_func .= "$perl $Bin/function/auto_eggNOG.pl -eggNOG -KEGG -GO $pep -workdir $function_each_dir \n";

        #$cmd_func .= "sh $function_dir/$name/STEP03_delete_tmp_files.sh\n";
        generateShell("$shell_dir/S5_function.$name.sh", $cmd_func ,\@shell_list);
    }
    close IN5;
}

if($step=~/6/){
    $Upload_dir = "$workdir/Upload";
    `mkdir -p $Upload_dir`;

    $pep_input = "$workdir/pep.list";
    open IN6, "<$pep_input" || die $!;

    while(<IN6>){
        chomp;
        my ($name,$pep)=split /\s+/,$_;
        $pep = abs_path($pep);

        my $Upload_each_dir = "$Upload_dir/$name";

        ###
        my $cmd_report = "mkdir -p $Upload_each_dir/01.Filter/FastQC \n";
        $cmd_report .= "ln -sf $filter_dir/$name/$name.raw.R1_fastqc/Images/per_base_quality.png $Upload_each_dir/01.Filter/FastQC/${name}_raw_1.png \n";
        $cmd_report .= "ln -sf $filter_dir/$name/$name.raw.R2_fastqc/Images/per_base_quality.png $Upload_each_dir/01.Filter/FastQC/${name}_raw_2.png \n";
        $cmd_report .= "ln -sf $filter_dir/$name/$name.qc_stat.txt $Upload_each_dir/01.Filter/${name}_qc_stat.txt \n";
        $cmd_report .= "$python $Bin/../../file/txt2excel.py $Upload_each_dir/01.Filter/${name}_qc_stat.txt \n";

        $cmd_report .= "ln -sf $filter_dir/$name/$name.clean.R1_fastqc/Images/per_base_quality.png $Upload_each_dir/01.Filter/FastQC/${name}_clean_1.png \n";
        $cmd_report .= "ln -sf $filter_dir/$name/$name.clean.R2_fastqc/Images/per_base_quality.png $Upload_each_dir/01.Filter/FastQC/${name}_clean_2.png \n";

        ###
        $cmd_report .= "mkdir -p $Upload_each_dir/02.Assembly \n";
        $cmd_report .= "ln -sf $assembly_dir/$name/$name.genome.fa $Upload_each_dir/02.Assembly/$name.genome.fa \n";
        $cmd_report .= "ln -sf $assembly_dir/$name/${name}_fa.stat.txt $Upload_each_dir/03.Assembly/${name}_fa.stat.txt \n";
        $cmd_report .= "python $Bin/../../file/txt2excel.py $Upload_each_dir/03.Assembly/${name}_filter_fa.stat.txt \n";

        # $cmd_report .= "ln -sf $assembly_each_dir/${name}_filter.length.png $Upload_each_dir/03.Assembly/${name}_filter.length.png \n";
        $cmd_report .= "ln -sf $assembly_dir/$name/$name.genome.fa.contigs.length_distribution.png $Upload_each_dir/02.Assembly/${name}_contigs_length_distribute.png \n";
        $cmd_report .= "ln -sf $assembly_dir/$name/$name.genome.fa.contigs.length_distribution.pdf $Upload_each_dir/02.Assembly/${name}_contigs_length_distribute.pdf \n";
        # $cmd_report .= "ln -sf $assembly_dir/$name/$name.contig.stat $Upload_each_dir/02.Assembly/${name}_contig_stat.txt \n";
        # $cmd_report .= "$python $Bin/../../file/txt2excel.py $Upload_each_dir/02.Assembly/${name}_contig_stat.txt \n";

        ###
        $cmd_report .= "mkdir -p $Upload_each_dir/03.Contig \n";
        # $cmd_report .= "ln -sf $contig_dir/$name/01.abundance/$name.sort.filter.cov.contig.abundance $Upload_each_dir/03.Contig/${name}_Contig_aubn.xls \n";
        # $cmd_report .= "ln -sf $contig_dir/$name/02.taxonomy/$name.contig.taxids.lca $Upload_each_dir/03.Contig/${name}_Contig_taxon.xls \n";
        $cmd_report .= "ln -sf $contig_dir/$name/${name}_Contig_abun_and_lineage.txt $Upload_each_dir/03.Contig/${name}_Contig_abun_and_lineage.txt \n";
        $cmd_report .= "$python $Bin/../../file/txt2excel.py $Upload_each_dir/03.Contig/${name}_Contig_abun_and_lineage.txt \n";

        ###
        $cmd_report .= "mkdir -p $Upload_each_dir/04.Gene \n";
        $cmd_report .= "ln -sf $gene_dir/$name/01.gene_prediction/$name/$name.gff  $Upload_each_dir/04.Gene/${name}.gff \n";
        $cmd_report .= "ln -sf $gene_dir/$name/01.gene_prediction/$name/$name.faa  $Upload_each_dir/04.Gene/${name}.faa \n";
        $cmd_report .= "ln -sf $gene_dir/$name/01.gene_prediction/$name/$name.ffn  $Upload_each_dir/04.Gene/${name}.ffn \n";

        $cmd_report .= "ln -sf $gene_dir/$name/01.gene_prediction/${name}.ffn.contigs.length_distribution.png $Upload_each_dir/04.Gene/${name}_gene_length_distribute.png \n";
        $cmd_report .= "ln -sf $gene_dir/$name/01.gene_prediction/${name}.ffn.contigs.length_distribution.pdf $Upload_each_dir/04.Gene/${name}_gene_length_distribute.pdf \n";
        $cmd_report .= "ln -sf $gene_dir/$name/02.gene_abun/$name.rpkm  $Upload_each_dir/04.Gene/${name}_gene_abun.txt \n";
        $cmd_report .= "$python $Bin/../../file/txt2excel.py  $Upload_each_dir/04.Gene/${name}_gene_abun.txt \n";


        ###
        $cmd_report .= "mkdir -p $Upload_each_dir/05.Function \n";
        $cmd_report .= "cd $Upload_each_dir/05.Function && mkdir -p CARD COG eggNOG GO KEGG VFDB \n";

        $cmd_report .= "ln -sf $function_dir/$name/CARD/card_anno.txt  $Upload_each_dir/05.Function/CARD/${name}_CARD.txt \n";
        $cmd_report .= "$python $Bin/../../file/txt2excel.py $Upload_each_dir/05.Function/CARD/${name}_CARD.txt \n";
        $cmd_report .= "ln -sf $function_dir/$name/VFDB/top1_VFDB_anno.txt  $Upload_each_dir/05.Function/VFDB/${name}_VFDB.txt \n";
        $cmd_report .= "$python $Bin/../../file/txt2excel.py $Upload_each_dir/05.Function/VFDB/${name}_VFDB.txt \n";
        
        $cmd_report .= "ln -sf $function_dir/$name/COG/all.COG.class.txt $Upload_each_dir/05.Function/COG/${name}_COG.txt \n";
        $cmd_report .= "$python $Bin/../../file/txt2excel.py $Upload_each_dir/05.Function/COG/${name}_COG.txt \n";

        $cmd_report .= "ln -sf $function_dir/$name/COG/all.COG.bar.png  $Upload_each_dir/05.Function/COG/${name}_COG.png \n";
        $cmd_report .= "ln -sf $function_dir/$name/eggNOG/eggNOG.txt $Upload_each_dir/05.Function/eggNOG/${name}_eggNOG.txt \n";
        $cmd_report .= "$python $Bin/../../file/txt2excel.py $Upload_each_dir/05.Function/eggNOG/${name}_eggNOG.txt \n";

        # go
        $cmd_report .= "if [ -e $function_dir/$name/GO_eggNOG/$name.faa.GO.png ];then \n";
        $cmd_report .= "\t ln -sf $function_dir/$name/GO_eggNOG/$name.faa.GO2Gene.txt $Upload_each_dir/05.Function/GO/${name}_GO2Gene.txt \n";
        $cmd_report .= "\t $python $Bin/../../file/txt2excel.py $Upload_each_dir/05.Function/GO/${name}_GO2Gene.txt \n";
        $cmd_report .= "\t ln -sf $function_dir/$name/GO_eggNOG/$name.faa.GO.png $Upload_each_dir/05.Function/GO/${name}_GO.png \n";
        $cmd_report .= "else\n\techo \"No GO anno\"\n";
        $cmd_report .= "fi\n";

        # KEGG
        $cmd_report .= "if [ -e $function_dir/$name/KEGG/$name.faa.KEGG.png ];then \n";
        $cmd_report .= "\t ln -sf $function_dir/$name/KEGG/$name.faa.kegg.txt  $Upload_each_dir/05.Function/KEGG/${name}_KEGG.txt \n";
        $cmd_report .= "\t $python $Bin/../../file/txt2excel.py $Upload_each_dir/05.Function/KEGG/${name}_KEGG.txt \n";
        $cmd_report .= "\t ln -sf $function_dir/$name/KEGG/$name.faa.KEGG.png  $Upload_each_dir/05.Function/KEGG/${name}_KEGG.png\n";
        $cmd_report .= "else\n\techo \"No KEGG anno\"\n";
        $cmd_report .= "fi\n";

        # $cmd_report .= "ln -sf $function_dir/$name/KEGG/$name.faa.kegg.txt  $Upload_each_dir/05.Function/KEGG/${name}_KEGG.txt \n";
        # $cmd_report .= "python $Bin/../../file/txt2excel.py $Upload_each_dir/05.Function/KEGG/${name}_KEGG.txt \n";
        # $cmd_report .= "ln -sf $function_dir/$name/KEGG/$name.faa.KEGG.png  $Upload_each_dir/05.Function/KEGG/${name}_KEGG.png \n";

        $cmd_report .= "$perl  $Bin/meta_denovo_report.pl -outdir $Upload_each_dir -name $name \n";

        if ($task_id){
            $cmd_report .= "mkdir -p $AnalysisResults/$task_id/$name \n";
            $cmd_report .= "cp -R $workdir/Upload/$name $AnalysisResults/$task_id/ \n";
            $cmd_report .= "cd $AnalysisResults/$task_id/ \n";
            $cmd_report .= "zip -qr $name.zip $name \n";
        }

        generateShell("$shell_dir/S6_report.$name.sh", $cmd_report ,\@shell_list);

    }
    close IN6;
}

open MAIN,">$shell_dir/main.sh" || die $!;
print MAIN join (" && \\\n",@shell_list) . "\n";
close MAIN;

`bash $shell_dir/main.sh`;

if ($task_id){
    my $cmd_log = "cp $shell_dir/*sh.e $shell_dir/*sh.o $log_dir \n";
    $cmd_log .= "cd $workdir \n";
    $cmd_log .= "zip -qr ${task_id}_log.zip log \n";
    system($cmd_log);
}

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
