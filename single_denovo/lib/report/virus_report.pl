use lib "/sdbb/share/lib/GDHR";
use strict;
use GDHR;
use Utils;


#定义报告路径
my $sample=shift;
my $andir=shift;
my $dir="$andir/$sample";

our @blank = map {"&nbsp;" x $_} 1 .. 20;
my $blank6 = blank(6);

#############################################################################
# 标题
my $report;
$report = GDHR->new(-outdir => $dir,
    -pipe                   => "单病毒有参拼接分析报告",
    -nonlazy                => 0);

#############################################################################
# 1.项目介绍
my $section;
$section = $report->section(id => "introduction");
$section->menu("项目概述");
$section->desc("$blank6 通过分离培养获得单病毒株，进行建库测序，获得较为纯净的毒株序列，本项目对单毒株（样本编号为$sample）的二代双端测序数据进行分析，通过输入参考基因组，进行比对拼接，得到一致性序列，并对该序列进行一系列生信分析。");

# 流程图
$section->menu("信息分析流程",-icon => "./src/image/icon/1.fenxi_liucheng.png");
$section->desc("$blank6 生物信息分析流程图，如下图所示:");
$section->img2html(
    -file => "./src/image/virus.png",
    -name => "信息分析流程图");

#############################################################################
# 2.数据质控
my $section;
$section = $report->section(id => "stat_show");
$section->menu("数据质控",-icon => "./src/image/icon/2.raw_qc.png");
FileLink($section, "$blank6  全部原始数据质控结果请点击  <link=./1.qc>");

#数据过滤介绍
$section->desc("$blank6 测序数据的产生是经过了 DNA 提取、建库、测序多个步骤的，这些步骤中产生的无效数据会对后续信息分析带来严重干扰，如测序接头序列，建库长度的偏差，以及测序错误、低质量碱基、未测出的碱基（以 N 表示）等情况，须通过一些手段将上述无效数据过滤掉，以保证分析的正常进行。");
$section->desc("$blank6 过滤标准主要有： 1) 过滤掉含有接头序列的 reads； 2) 截掉 reads 部分区域的低质量碱基； 3）去除 adapter 序列 ");

# 2.1过滤前
$section->submenu("测序数据质量分析",-icon => "./src/image/icon/2.1guolv.png");
$section->ssubmenu("过滤前碱基质量图");
#双端
if (-s  "${dir}/1.qc/${sample}_1_fastqc/Images/per_base_quality.png" ){
$section->desc("$blank6 样本 read_1和read2端过滤前，碱基质量图如下：");
my @name = ["imgs2html2"];
$section->img2html2(
    -file1 => "./1.qc/${sample}_1_fastqc/Images/per_base_quality.png",
    -name1  => "read_1端过滤前碱基质量图",
    -file2 => "./1.qc/${sample}_2_fastqc/Images/per_base_quality.png",
    -name2  => "read_2端过滤前碱基质量图",
    -names  => "img");


FileLink($section, "$blank6  更多read_1 过滤前数据质控报告请点击  <link=./1.qc/${sample}_1_fastqc.html>");
FileLink($section, "$blank6  更多read_2 过滤前数据质控报告请点击  <link=./1.qc/${sample}_2_fastqc.html>");

#单端
}elsif (-s "${dir}/1.qc/${sample}_fastqc/Images/per_base_quality.png" ){
    $section->desc("$blank6 样本${sample}过滤前碱基质量图如下：");
    $section->img2html(
        -file => "./1.qc/${sample}_fastqc/Images/per_base_quality.png",
        -name  => "过滤前碱基质量图");
    FileLink($section, "$blank6  更多过滤前数据质控报告请点击  <link=./1.qc/${sample}_fastqc.html>");
}

#质控图说明
$section->desc("$blank6 说明：横坐标表示测序位置，纵坐标为测序质量值图中，横轴代表位置，纵轴 quality。红色表示中 位数，黄色是 25%-75%区间，触须是 10%-90%区间，蓝线是平均数。随着测序的进行，酶的活性会逐步下降，因此到达一定测序长度后，碱基质量值也会随之下降。");


# 2.2过滤后
$section->ssubmenu("过滤后碱基质量图");
#双端
if (-s  "${dir}/1.qc/${sample}_1.clean_fastqc/Images/per_base_quality.png" ){
$section->desc("$blank6 样本 read_1和read2端过滤后，碱基质量图如下：");
my @name = ["imgs2html2"];
$section->img2html2(
    -file1 => "./1.qc/${sample}_1.clean_fastqc/Images/per_base_quality.png",
    -name1  => "read_1端过滤后碱基质量图",
    -file2 => "./1.qc/${sample}_2.clean_fastqc/Images/per_base_quality.png",
    -name2  => "read_2端过滤后碱基质量图",
    -names  => "img");


FileLink($section, "$blank6  更多read_1 过滤后数据质控报告请点击  <link=./1.qc/${sample}_1.clean_fastqc.html>");
FileLink($section, "$blank6  更多read_2 过滤后数据质控报告请点击  <link=./1.qc/${sample}_2.clean_fastqc.html>");

#单端
}elsif (-s "${dir}/1.qc/${sample}.clean_fastqc/Images/per_base_quality.png" ){
    $section->desc("$blank6 样本${sample}过滤后碱基质量图如下：");
    $section->img2html(
        -file => "./1.qc/${sample}.clean_fastqc/Images/per_base_quality.png",
        -name  => "过滤后碱基质量图");
    FileLink($section, "$blank6  更多过滤后数据质控报告请点击  <link=./1.qc/${sample}.clean_fastqc.html>");
}


#质控图说明
$section->desc("$blank6 说明：横坐标表示测序位置，纵坐标为测序质量值图中，横轴代表位置，纵轴 quality。红色表示中 位数，黄色是 25%-75%区间，触须是 10%-90%区间，蓝线是平均数。随着测序的进行，酶的活性会逐步下降，因此到达一定测序长度后，碱基质量值也会随之下降。");



#############################################################################
#1.2数据统计
$section->submenu("数据统计",-icon => "./src/image/icon/2.2stat.png");

#常规指标
if (-s  "${dir}/1.qc/${sample}.basic.stat.txt" ){
$section->tsv2html(
    -file   => "${dir}/1.qc/${sample}.basic.stat.txt",
    -name   => "数据统计结果",
    -header => 1);
}else{$section->desc("$blank6 无数据统计结果");}





############################################################################
#2有参拼接
if (-s "$dir/2.ass/${sample}_consensus.fa"){
    $section = $report->section(id => "ass");
    $section->menu("有参拼接",-icon => "./src/image/icon/3.ass.png");
    FileLink($section, "$blank6 全部有参拼接结果请点击  <link=./2.ass>");

    #拼接contigs序列
    $section->submenu("拼接序列",-icon => "./src/image/icon/3.1dna.png");
    FileLink($section, "$blank6 拼接序列请点击  <link=./2.ass/${sample}_consensus.fa>");


    $section->submenu("拼接序列质控",-icon => "./src/image/icon/3.2ass_qc.png");

    #深度覆盖度评估
    $section->ssubmenu("拼接基因组覆盖度及深度统计");
    $section->desc("$blank6 下图展示原始数据比对上拼接序列的覆盖度，测序深度（5x，30x，100x），和一致性");

    if (-s "$dir/2.ass/bam.stat.txt"){
    $section->tsv2html(
        -file   => "$dir/2.ass/bam.stat.txt",
        -name   => "拼接基因组覆盖度与测序深度统计",
        -header => 1);
    FileLink($section, "$blank6  拼接结果详情请点击  <link=./2.ass>");

    }else{
        $section->desc("$blank6 无拼接基因组覆盖度与测序深度统计结果");
    }


    $section->ssubmenu("组装基因深度分布图");
    $section->img2html(
        -file => "./2.ass/bed.depth.stat.png",
        -name => "组装基因深度分布图");
    $section->desc("注： 散点图，横坐标为GC含量，纵坐标为reads覆盖深度；两侧散点图，分别为GC含量、测序深度的滑窗频数分布。通过GC_depth分布图可以看出测序是否有明显的GC偏向，也可以判断样品是否存在污染等情况。如果多数的点集中分布在一个比较窄的范围内，则表明样本不存在污染；如果分布在多个区域，则表明样本中可能存在其它物种的污染。当存在污染时，可根据点分布的集中区判断污染程度，例如根据GC分布大致判断污染物种的数量，或者根据污染序列的测序reads覆盖深度来大致推测污染的比例等。");


    #################################################################################################
    #4变异注释
    $section = $report->section(id => "variants");
    $section->menu("变异统计",-icon => "./src/image/icon/4.3shijunti.png");

    if (-s "$dir/2.ass/variants.txt"){
    $section->tsv2html(
        -file   => "$dir/2.ass/variants.txt",
        -name   => "变异信息统计",
        -header => 1);

    $section->desc("$blank6 说明：REF_NAME：参考序列名；POS：变异位点；REF:参考序列对应位置的碱基；ALT：拼接序列对应位置变异后的碱基；TOTAL_DP：对应位置测序深度 ");
    FileLink($section, "$blank6  拼接结果详情请点击  <link=./2.ass/variants.xlsx>");
    }else{
        $section->desc("$blank6 无变异信息统计结果");
    }



# 无参组装
}elsif (-e "${dir}/2.ass/${sample}_scaffolds.fasta"){
$section = $report->section(id => "assemble");
$section->menu("组装基因组",-icon => "./src/image/icon/3.ass.png");

    if (-e "${dir}/2.ass/${sample}_scaffolds.fasta"){
        #组装序列
        FileLink($section, "$blank6  全部组装质控结果请点击  <link=./2.ass>");
        $section->submenu("组装序列",-icon => "./src/image/icon/3.1dna.png");
        FileLink($section, "$blank6  过滤前组装scaffolds序列请点击  <link=./2.ass/scaffolds.fasta>");
        FileLink($section, "$blank6  过滤后组装scaffolds序列请点击  <link=./2.ass/${sample}_scaffolds.fasta>");



        $section->submenu("组装质控",-icon => "./src/image/icon/3.2ass_qc.png");
        #####################
        ### 3.2.1 基本组装指标
        $section->ssubmenu("基本组装指标");
        $section->desc("$blank6 组装完成后，可得到一系列用于表征组装质量的参数，包括scaffolds数目，scaffolds总长度、N50等。一般来说，scaffolds数目越少，N50、N90值越高，表示组装效果越好。");
        if (-e "${dir}/2.ass/${sample}_scaffolds_fa.stat.txt"){
        $section->tsv2html(
            -file   => "${dir}/2.ass/${sample}_scaffolds_fa.stat.txt",
            -name   => "${sample}组装质控信息",
            -header => 1);
        }else{
            $section->desc("$blank6 无组装质控信息");}
        $section->desc("$blank6 注：N50：将组装序列按照长度，从长到短进行排序，将序列长度依次相加，当累计长度达到组装序列总长度50%时，当前累加的序列长度即为 N50，累加的contig的个数即为L50，同理可得N75与N90，L75与L90。");



        #######################
        # 3.2.2组装序列长度分布图
        $section->ssubmenu("组装基因组长度分布图");
        if (-e "${dir}/2.ass/${sample}_scaffolds.length.png"){
        $section->img2html(
            -file => "./2.ass/${sample}_scaffolds.length.png",
            -name => "组装序列长度分布图");
        }else{
            $section->desc("$blank6 无组装序列长度分布图");}



        #############
        # 3.2.3一致性
        $section->ssubmenu("组装基因组一致性");
        $section->desc("$blank6 将抽取reads比对回组装序列上，统计抽取reads对序列的覆盖情况，从而评估组装结果与测序数据的一致性。较高的覆盖率（95%以上）认为组装结果和reads有比较好的一致性。");
        if (-s  "${dir}/2.ass/uniformity.txt"){
            $section->tsv2html(
                -file   => "${dir}/2.ass/uniformity.txt",
                -name   => "组装一致性统计结果",
                -header => 1);
            FileLink($section, "$blank6  每条scaffolds序列的组装一致性统计结果请点击  <link=./2.ass/uniformity_all.txt>");
        }else{
            $section->desc("$blank6 无组装一致性统计结果");}



        #####################
        # 3.2.4插入片段分布图
        $section->ssubmenu("插入片段分布图");
        if (-e "${dir}/2.ass/insertsize.png"){
            $section->img2html(
                -file => "./2.ass/insertsize.png",
                -name => "insert size分布图");
        }else{
            $section->desc("$blank6 无插入片段分布图或原始序列为单端");
        }



        #####################
        # 3.2.5 GC分布
        $section->ssubmenu("组装基因深度分布");
        $section->desc("$blank6 通过用滑窗法切割组装基因组，来评估测序是否有明显的GC偏向，也可判断样品是否存在污染等情况。如果多数的点集中分布在一个比较窄的范围内，则表明样本不存在污染；如果分布在多个区域，则表明样本中可能存在其它物种的污染。当存在污染时，可根据点分布的集中区判断污染程度，例如根据GC分布大致判断污染物种的数量，或者根据污染序列的测序reads覆盖深度来大致推测污染的比例等");
        $section->img2html(
            -file => "./2.ass/depth_base.stat.depth_GC.png",
            -name => "组装基因深度分布");
        $section->desc("注： 散点图，横坐标为GC含量，纵坐标为reads覆盖深度；两侧散点图，分别为GC含量、测序深度的滑窗频数分布。通过GC_depth分布图可以看出测序是否有明显的GC偏向。");



        #####################
        # 3.2.6 checkm
        # $section->ssubmenu("组装基因组完整度/污染度评估");
        # $section->desc("$blank6 一个谱系(lineage)的基因组会有一些共享的基因，称之为lineage-specific marker set。");
        # $section->desc("$blank6 评估一个基因组的完整度，首先判断这个基因组的lineage，通过检索该lineage的marker set，如果marker set的基因全部找到了，则认为这个基因组完整度100%，若检测到多个marker，则有可能是杂合或者污染。");
        # $section->desc("$blank6 本次分析对单菌基因组组装的完整度、杂合度质量评估，以及根据Marker基因鉴定基因组的系统分类");

        # if (-s  "${dir}/2.ass/check.txt"){
        #     $section->tsv2html(
        #         -file   => "${dir}/2.ass/check.txt",
        #         -name   => "组装完整度和准确度统计",
        #         -header => 1);
        # }else{
        #     $section->desc("$blank6 无组装完整度和准确度统计结果");}



        #####################
        # 3.2.7比对参考基因组
        if (-s  "${dir}/2.ass/ref_cov_stat.txt"){
            $section->ssubmenu("组装基因组比对参考基因组");
            $section->tsv2html(
                -file   => "${dir}/2.ass/ref_cov_stat.txt",
                -name   => "比对参考基因组统计",
                -header => 1);
        }




    }else{
        $section->desc("本次分析组装失败，请检查序列信息。");}



}





#参考文献
$section = $report->section(id => "pubmed");
$section->menu("参考文献",-icon => "./src/image/icon/6.pub.png");
my $desc = <<REF;
<p>
<ol>
<li>Alcock, B.P. et al. CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. Nucleic Acids Res 48, D517-D525 (2020).</li>
<li>Huerta-Cepas, J. et al. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Res 47, D309-D314 (2019).</li>
<li>Langmead, B., Trapnell, C., Pop, M. & Salzberg, S.L. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology 10, R25 (2009).</li>
<li>Leggett, R.M., Ramirez-Gonzalez, R.H., Clavijo, B.J., Waite, D. & Davey, R.P. Sequencing quality assessment tools to enable data-driven informatics for high throughput genomics. Front Genet 4, 288-288 (2013).</li>
<li>Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34, 3094-3100 (2018).</li>
<li>Liu, B., Zheng, D., Jin, Q., Chen, L. & Yang, J. VFDB 2019: a comparative pathogenomic platform with an interactive web interface. Nucleic Acids Res 47, D687-D692 (2019).</li>
<li>Parks, D.H., Imelfort, M., Skennerton, C.T., Hugenholtz, P. & Tyson, G.W. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research 25, 1043-1055 (2015).</li>
<li>alzberg, S.L. et al. GAGE: A critical evaluation of genome assemblies and assembly algorithms. Genome research 22, 557-567 (2012).</li>
<li>Seemann, T. Prokka: rapid prokaryotic genome annotation. Bioinformatics 30, 2068-2069 (2014).
<li>Simão, F.A., Waterhouse, R.M., Ioannidis, P., Kriventseva, E.V. & Zdobnov, E.M. BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics 31, 3210-3212 (2015).</li>
<li>Forouzan, E., Shariati, P., Mousavi Maleki, M.S., Karkhane, A.A. & Yakhchali, B. Practical evaluation of 11 de novo assemblers in metagenome assembly. Journal of Microbiological Methods 151, 99-105 (2018).</li>
</ol>
</p>
REF

$section->add_html($desc);


$report->write();

