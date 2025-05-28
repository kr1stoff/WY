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
    -pipe                   => "单真菌全基因组组装分析报告",
    -nonlazy                => 0);

#############################################################################
# 1.项目介绍
my $section;
$section = $report->section(id => "introduction");
$section->menu("项目概述");
$section->desc("$blank6 通过分离培养获得单真菌株，进行建库测序，获得较为纯净的菌株序列，本项目对单细菌株（样本编号为$sample）的二代双端测序数据进行分析，通过denovo拼接，最终得到高质量的基因组序列，并对该序列进行基因预测及基因注释。");

# 流程图
$section->menu("信息分析流程",-icon => "./src/image/icon/1.fenxi_liucheng.png");
$section->desc("$blank6 生物信息分析流程图，如下图所示:");
$section->img2html(
    -file => "./src/image/fungi.png",
    -name => "分析流程图");

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
#2.2测序质量评估及质控
$section->submenu("数据统计",-icon => "./src/image/icon/2.2stat.png");

#常规指标
if (-s  "${dir}/1.qc/${sample}.basic.stat.txt" ){
$section->tsv2html(
    -file   => "${dir}/1.qc/${sample}.basic.stat.txt",
    -name   => "数据统计结果",
    -header => 1);
}else{$section->desc("$blank6 无数据统计结果");}





#############################################################################
#3序列拼接
$section = $report->section(id => "ass");
$section->menu("组装基因组",-icon => "./src/image/icon/3.ass.png");
FileLink($section, "$blank6  全部组装质控结果请点击  <link=./2.ass>");



if (-e "${dir}/2.ass/${sample}_scaffolds.fasta"){

    #组装序列
    $section->submenu("组装序列",-icon => "./src/image/icon/3.1dna.png");
    FileLink($section, "$blank6  过滤前组装scaffolds序列请点击  <link=./2.ass/scaffolds.fasta>");
    FileLink($section, "$blank6  过滤后组装scaffolds序列请点击  <link=./2.ass/${sample}_scaffolds.fasta>");

    $section->submenu("组装质控",-icon => "./src/image/icon/3.2ass_qc.png");
    $section->ssubmenu("基本组装指标");
    $section->desc("$blank6 组装完成后，可得到一系列用于表征组装质量的参数，包括scaffolds数目，scaffolds总长度、N50等。一般来说，scaffolds数目越少，N50、N90值越高，表示组装效果越好。");
    if (-e "${dir}/2.ass/${sample}_scaffolds_fa.stat.txt"){
    $section->tsv2html(
        -file   => "${dir}/2.ass/${sample}_scaffolds_fa.stat.txt",
        -name   => "${sample}组装质控信息",
        -header => 1);
    }else{
        $section->desc("$blank6 无组装质控信息");
    }

    $section->desc("$blank6 注：N50：将组装序列按照长度，从长到短进行排序，将序列长度依次相加，当累计长度达到组装序列总长度50%时，当前累加的序列长度即为 N50，同理可得N75与N90。");

    #组装基因组长度分布图
    $section->ssubmenu("组装基因组长度分布图");
    if (-e "${dir}/2.ass/${sample}_scaffolds.length.png"){
    $section->img2html(
        -file => "./2.ass/${sample}_scaffolds.length.png",
        -name => "组装基因组长度分布图");
    }else{
        $section->desc("$blank6 无组装基因组长度分布图");
    }


    #一致性
    $section->ssubmenu("组装基因组一致性");
    $section->desc("$blank6 将抽取reads比对回组装序列上，统计抽取reads对序列的覆盖情况，从而评估组装结果与测序数据的一致性。较高的覆盖率（95%以上）认为组装结果和reads有比较好的一致性。");
    if (-s  "${dir}/2.ass/uniformity.txt"){
    $section->tsv2html(
        -file   => "${dir}/2.ass/uniformity.txt",
        -name   => "组装基因组一致性统计结果",
        -header => 1);
    FileLink($section, "$blank6  每条scaffolds序列的组装一致性统计结果请点击  <link=./2.ass/uniformity_all.txt>");
    }else{
        $section->desc("$blank6 无组装基因组一致性统计结果");
    }

    #插入片段分布图
    $section->ssubmenu("插入片段分布图");
    if (-e "${dir}/2.ass/insertsize.png"){
    $section->img2html(
        -file => "./2.ass/insertsize.png",
        -name => "insert size分布图");
    }else{
        $section->desc("$blank6 无插入片段分布图或原始序列为单端");
    }

    #GC分布
    $section->ssubmenu("组装基因深度分布");
    $section->desc("$blank6 通过用滑窗法切割组装基因组，来评估测序是否有明显的GC偏向，也可判断样品是否存在污染等情况。如果多数的点集中分布在一个比较窄的范围内，则表明样本不存在污染；如果分布在多个区域，则表明样本中可能存在其它物种的污染。当存在污染时，可根据点分布的集中区判断污染程度，例如根据GC分布大致判断污染物种的数量，或者根据污染序列的测序reads覆盖深度来大致推测污染的比例等");
    $section->img2html(
        -file => "./2.ass/depth_base.stat.depth_GC.png",
        -name => "组装基因深度分布");
    $section->desc("注： 散点图，横坐标为GC含量，纵坐标为reads覆盖深度；两侧散点图，分别为GC含量、测序深度的滑窗频数分布。通过GC_depth分布图可以看出测序是否有明显的GC偏向。");

    #checkm
    $section->ssubmenu("组装完整度/污染度评估");
    $section->desc("$blank6 一个谱系(lineage)的基因组会有一些共享的基因，称之为lineage-specific marker set。");
    $section->desc("$blank6 评估一个基因组的完整度，首先判断这个基因组的lineage，通过检索该lineage的marker set，如果marker set的基因全部找到了，则认为这个基因组完整度100%，若检测到多个marker，则有可能是杂合或者污染。");
    $section->desc("$blank6 本次分析对单菌基因组组装的完整度、杂合度质量评估，以及根据Marker基因鉴定基因组的系统分类");

    if (-s  "${dir}/2.ass/check.txt"){
    $section->tsv2html(
        -file   => "${dir}/2.ass/check.txt",
        -name   => "组装完整度和准确度统计",
        -header => 1);
    }else{
        $section->desc("$blank6 无组装完整度和准确度统计结果");
    }

    # 比对参考基因组
    if (-s  "${dir}/2.ass/ref_cov_stat.txt"){
    $section->ssubmenu("组装基因组比对参考基因组");
    $section->tsv2html(
        -file   => "${dir}/2.ass/ref_cov_stat.txt",
        -name   => "比对参考基因组统计",
        -header => 1);
    }
}




    ############################################################################################################################################
    #4基因预测
    $section = $report->section(id => "predict");
    $section->menu("基因预测",-icon => "./src/image/icon/4.pre.png");
    FileLink($section, "$blank6 全部基因预测结果请点击  <link=./3.pre>");


    if (-e "$andir/3.pre/predict/genemark.gtf"){

        # 基因长度分布图
        $section->submenu("预测基因长度分布",-icon => "./src/image/icon/4.1length.png");
        $section->img2html(
            -file => "./3.pre/predict.length.png",
            -name => "预测基因的长度分布图");
        $section->desc("$blank6 横坐标为基因长度范围（bp），纵坐标为对应长度范围的预测的基因数量");


        # 预测基因类型统计
        $section->submenu("预测基因分类统计",-icon => "./src/image/icon/4.2class.png");
        if (-s "${dir}/3.pre/predict_kind.txt"){
        $section->tsv2html(
            -file   => "${dir}/3.pre/predict_kind.txt",
            -name   => "预测基因分类统计表",
            -header => 1);
        }else{
            $section->desc("$blank6 无预测基因分类统计表结果");
        }

        # 4.3.3重复序列预测
        $section->submenu("重复序列预测",-icon => "./src/image/icon/4.5repeat.png");
        $section->desc("$blank6 生物细胞中的DNA序列里面包含许多重复序列（repeated sequence），主要可分为两大类，分别是串联重复序列（也叫串接重复序列，Tandem repeat）与散置重复序列（Interspersed repeat）。 串联重复和散置重复的区别在于重复的部分是否相邻分布，相邻就是串联重复，不相邻就是散置重复。");

        $section->ssubmenu("串联重复序列预测");
        if (-s "${dir}/3.pre/island.txt"){
            $section->tsv2html(
                -file   => "${dir}/3.pre/TandemRepeat.txt",
                -name   => "串联重复序列预测结果",
                -header => 1);
            FileLink($section, "$blank6  全部串联重复序列预测结果请点击  <link=./3.pre/TandemRepeat.txt>");
        }else{
            $section->desc("$blank6 无串联重复序列预测结果");
        }

        section->ssubmenu("散置重复序列预测");
        if (-s "${dir}/3.pre/INEs.txt"){
            $section->tsv2html(
                -file   => "${dir}/3.pre/INEs.txt",
                -name   => "散置重复序列预测结果",
                -header => 1);
            FileLink($section, "$blank6  全部散置重复序列预测结果请点击  <link=./3.pre/INEs.txt>");
        }else{
            $section->desc("$blank6 无散置重复序列预测结果");
        }




        ################################################################################################################
        #5基因注释
        $section = $report->section(id => "anno");
        $section->menu("基因注释",-icon => "./src/image/icon/5.anno.png");
        FileLink($section, "$blank6  全部基因注释结果请点击  <link=./4.anno>");

        #5.1 毒力因子注释
        $section->submenu("毒力因子注释",-icon => "./src/image/icon/5.1vir.png");
        $section->desc("$blank6 本分析使用通用的毒力数据库，用于专门研究致病细菌、衣原体和支原体致病因子的数据库。其包含 75 个属，共 1800 个致病因子，30053 个与毒力因子相关的基因。");
        $section->desc("$blank6 本项目将预测的氨基酸序列，与上述数据库进行比对，把预测基因和其相对应的毒力因子功能注释信息结合起来，得到注释结果。");
        $section->desc("$blank6 毒力因子包括细菌毒素,调节细菌附着的细胞表面蛋白,保护细菌细胞表面的碳水化合物和蛋白质和水 解酶等。在一些物种中，不同的独立因子组合能够引起不同的疾病。进一步的比较基因组分析能够发现更 多关于生物和病原体进化的信息。");

        # 毒力基因详细结果
        if (-s "${dir}/4.anno/Virulence.txt"){
            $section->tsv2html(
                -file   => "${dir}/4.anno/Virulence.txt",
                -name   => "毒力基因比对结果",
                -header => 1);
            FileLink($section, "$blank6  全部毒力基因注释结果点击  <link=./4.anno/Virulence.txt>");
            FileLink($section, "$blank6  毒力基因注释序列点击  <link=./4.anno/Virulence_gene.fa>");
            FileLink($section, "$blank6  毒力基因参考序列点击  <link=./4.anno/Virulence_ref.fa>");
        }else{
            $section->desc("$blank6 无毒力基因比对结果");}


        ################
        # 5.2 耐药基因注释 
        $section->submenu("耐药基因注释",-icon => "./src/image/icon/5.2drug.png");
        
        $section->ssubmenu("CARD注释");
        $section->desc("$blank6 耐药基因注释选用通用的数据库进行比对，有 2359 条序列，3567 条 Antibiotic Resitance Ontology Term。");
        $section->desc("$blank6 将预测的的氨基酸序列，与上述数据库进行比对，把预测的基因和其相对应的耐药功能注释信息结合起来，得到注释结果。");
        if (-s "${dir}/4.anno/CARD_anno.txt"){
            $section->tsv2html(
                -file   => "${dir}/4.anno/CARD_anno.txt",
                -name   => "耐药基因注释",
                -header => 1);
            FileLink($section, "$blank6  耐药基因注释结果点击  <link=./4.anno/CARD_anno.txt>");
            FileLink($section, "$blank6  耐药基因注释序列点击  <link=./4.anno/CARD_contig.fa>");
            FileLink($section, "$blank6  耐药基因注释参考序列点击  <link=./4.anno/CARD_ref.fa>");
        }else{
            $section->desc("$blank6 无CARD注释结果");
        }
        ####
        $section->ssubmenu("Resfinder注释");
        $section->desc("$blank6 耐药基因注释选用通用的数据库进行比对，有 2359 条序列，3567 条 Antibiotic Resitance Ontology Term。");
        $section->desc("$blank6 将预测的的氨基酸序列，与上述数据库进行比对，把预测的基因和其相对应的耐药功能注释信息结合起来，得到注释结果。");
        if (-s "${dir}/4.anno/Resfinder_anno.txt"){
            $section->tsv2html(
                -file   => "${dir}/4.anno/Resfinder_anno.txt",
                -name   => "耐药基因注释",
                -header => 1);
            FileLink($section, "$blank6  耐药基因注释结果点击  <link=./4.anno/Resfinder_anno.txt>");
            FileLink($section, "$blank6  耐药基因注释序列点击  <link=./4.anno/ResFinder_Hit_in_genome_seq.fsa>");
        }else{
            $section->desc("$blank6 无Resfinder注释结果");


        # CAZy 注释
        $section->submenu("CAZy注释",-icon => "./src/image/icon/5.3cazy.png");
        $section->desc("$blank6 碳水化合物亦称糖类化合物，是自然界存在最多、分布最广的一类重要有机化合物，是一切生物体维持生命活动所需能量的主要来源。作用于各种糖复合物、寡糖和多糖等碳水化合物的酶类构成了地球上结构最多样的蛋白质集合。");
        $section->desc("$blank6 CAZY 全称为Carbohydrate-Active enZYmes Database，是碳水化合物酶相关的专业数据库，包括催化碳水化合物和糖复合物的生物合成、降解以及修饰的相关酶系家族。");
        $section->desc("$blank6 根据蛋白质结构域中氨基酸序列的相似性，分成六大类家族:糖苷水解酶类(GHs),糖苷转移酶类(GTs),多糖裂解酶类(PLs),糖水化合物酯酶类(CEs),辅助模块酶类(AAs),碳水化合物结合模块(CBMs)。");
        if (-s "${dir}/4.anno/CAZY.txt"){
        $section->img2html(
            -file => "./4.anno/CAZy_anno_stats.png",
            -name => "CAZy注释分类图");
        $section->tsv2html(
            -file   => "${dir}/4.anno/CAZY.txt",
            -name   => "CAZy注释",
            -header => 1);
        FileLink($section, "$blank6 全部CAZy注释结果点击  <link=./4.anno/CAZY.txt>");
        }else{
            $section->desc("$blank6 无CAZy注释结果");
        }



        # COG注释
        $section->submenu("COG注释",-icon => "./src/image/icon/5.4cog.png");
        $section->desc("$blank6 COG的中文释义即“同源蛋白簇”,COG（Clusters of Orthologous Groups ）注释是差异基因功能注释的一种方法。COG是由NCBI创建并维护的蛋白数据库，是对基因产物进行同源分类，较早的识别直系同源基因的数据库，通过对多种生物的蛋白质序列大量比较而来。通过比对可以将某个蛋白序列注释到某一个COG中，每一簇COG由直系同源序列构成，从而可以推测该序列的功能。COG数据库按照功能一共可以分为二十六类。");
        $section->desc("$blank6 COG注释作用：1. 通过已知蛋白对未知序列进行功能注释； 2. 通过查看指定的COG编号对应的protein数目，存在及缺失，从而能推导特定的代谢途径是否存在； 3. 每个COG编号是一类蛋白，将query序列和比对上的COG编号的proteins进行多序列比对，能确定保守位点，分析其进化关系");
        #COG注释分类图
        $section->img2html(
            -file => "./4.anno/all.COG.bar.png",
            -name => "COG注释分类图");


        $section->desc("说明：<br>A：RNA加工修饰;B：染色质结构和动力学;C：能量生成和转换;D：细胞周期控制、细胞分裂和染色体分裂;E：氨基酸转运代谢;F：核苷酸转运和代谢;G：碳水化合物转运代谢;H：辅酶转运和代谢;I: 脂肪转运代谢;J：翻译，核糖体结构和生物合成;K：转录;L：复制，重组和修复;M：细胞壁/膜/被膜的生物合成;N：细胞运动;O：翻译后修饰，蛋白质折叠和伴侣蛋白;P: 无机离子转运代谢;Q：次级代谢物生物合成，转运和代谢;R：主要功能预测;S：未知功能;T：信号转导机制;U： 胞内转运、分泌和小泡运输;V：抵御机制;W：胞外结构;X：动员组：噬菌体原，转座子;Y：核酸结构;Z：细胞骨架");


        # GO 注释
        $section->submenu("GO注释",-icon => "./src/image/icon/5.5go.png");
        $section->desc("$blank6 GO的全称是基因本体论(Gene Ontology)，提供了对基因功能与基因产物最为全面的描述。GO term按照分类，一共三个本体(ontology)，分别描述基因的分子功能(molecular function，MF)、细胞组分(cellular component，CC)、参与的生物过程(biological process，BP)。对于一个基因或蛋白来说，最终一个基因的产物会被多个GO术语进行注释。");
        $section->desc("$blank6 ");
        #GO注释分类图
        if (-s "${dir}/4.anno/GO_anno_stats.xls"){
        $section->tsv2html(
            -file   => "${dir}/4.anno/GO_anno_stats.xls",
            -name   => "GO分类注释结果",
            -header => 1);
        $section->img2html(
            -file => "./4.anno/GO_anno_stats_level2.png",
            -name => "GO注释分类图");
        $section->desc("说明：细胞组成(cellular component，CC)：一般用来描述基因产物的发挥作用的位置；生物过程(biological process，BP)：描述的是指基因产物所联系的一个大的生物功能，或者说是它们要完成的一个大的生物目标；分子功能(Molecular Function，MF)：主要是指基因产物分子所执行的任务。");
        FileLink($section, "$blank6  详细GO注释结果点击  <link=./4.anno/GO_anno_stats.xls>");
        }else{
            $section->desc("$blank6 无GO注释结果");
        }


        # KEGG 注释
        $section->submenu("KEGG注释",-icon => "./src/image/icon/5.6kegg.png");
        $section->desc("$blank6 基因会参与人体的各个通路。KEGG就是基于人体通路而形成的其中一个数据库，由日本京都大学生物信息学中心的Kanehisa实验室于1995年建立,数据库大致分为系统信息、基因组信息和化学信息三大类。");
        $section->desc("$blank6 ");
        # GO注释分类图
        if (-s "${dir}/4.anno/KEGG_anno.stat.txt"){
        $section->tsv2html(
            -file   => "${dir}/4.anno/KEGG_anno.txt",
            -name   => "KEGG注释结果",
            -header => 1);
        $section->img2html(
            -file => "./4.anno/KEGG_anno_stats.png",
            -name => "KEGG注释分类图");
        FileLink($section, "$blank6  详细KEGG注释分类结果点击  <link=./4.anno/KEGG_anno.stat.txt>");
        }else{
            $section->desc("$blank6 无KEGG注释结果");
        }

        # swiss-prot 注释
        $section->submenu("swiss-prot注释",-icon => "./src/image/icon/5.7swiss.png");
        $section->desc("$blank6 SWISS-PROT是经过注释的蛋白质序列数据库，由欧洲生物信息学研究所(EBI)维护。数据库由蛋白质序列条目构成，每个条目包含蛋白质序列、引用文献信息、分类学信息、注释等，注释中包括蛋白质的功能、转录后修饰、特殊位点和区域、二级结构、四级结构、与其它序列的相似性、序列残缺与疾病的关系、序列变异体和冲突等信息。");
        if (-s "${dir}/4.anno/swissprot_result.tsv"){
        $section->tsv2html(
            -file   => "${dir}/4.anno/swissprot_result.tsv",
            -name   => "swiss-prot注释",
            -header => 1);
        FileLink($section, "$blank6 全部swiss-prot注释结果点击  <link=./4.anno/swissprot_result.tsv>");

        }else{
            $section->desc("$blank6 无swiss-prot注释结果");
        }

        # # pfam 注释
        # $section->submenu("pfam注释");
        # $section->desc("$blank6 Pfam是蛋白质家族的数据库，根据多序列比对结果和隐马尔可夫模型，将蛋白质分为不同的家族。");
        # if (-s "${dir}/4.anno/pfam.txt"){
        # $section->tsv2html(
        #     -file   => "${dir}/4.anno/pfam.txt",
        #     -name   => "pfam注释",
        #     -header => 1);
        # FileLink($section, "$blank6 全部pfam注释结果点击  <link=./4.anno/pfam.txt>");
        # }else{
        #     $section->desc("$blank6 无pfam注释结果");
        # }



        ###########
        # 质粒注释
        $section->submenu("质粒注释",-icon => "./src/image/icon/5.8plasmid.png");
        $section->desc("$blank6 ");
        if (-s "${dir}/4.anno/plasmid.txt"){
            $section->tsv2html(
                -file   => "${dir}/4.anno/plasmid.txt",
                -name   => "质粒注释",
                -header => 1);
            FileLink($section, "$blank6 全部质粒注释结果点击  <link=./4.anno/plasmid.txt>");
        }else{
            $section->desc("$blank6 无质粒注释结果");}



    }else{
        $section->desc("本次分析预测失败，请检查序列信息。");
    }


}else{
    $section->desc("本次分析组装失败，请检查序列信息。");
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

