use lib "/sdbb/share/lib/GDHR";
use strict;
use GDHR;
use Utils;

# 读入命令行参数
my $sample=shift;
my $dir=shift;
our @blank = map {"&nbsp;" x $_} 1 .. 20;
my $blank6 = blank(6);

# 标题
my $report;
$report = GDHR->new(-outdir => $dir,
    -pipe                   => "溯源分型与进化树分析报告",
    -nonlazy                => 1);

#进化树
my $section;
$section = $report->section(id => "introduction");
$section->menu("进化树介绍");
$section->desc("$blank6 进化树，又称系统发育树，展示具有亲缘关系的物种/基因之间的种系发生历史的树状图。进化树由结点和进化分支组成，每一结点表示一个分类学单元，进化分支定义了分类单元之间的关系。");

# 结果
my $section;
$section = $report->section(id => "result");
$section->menu("进化树结果");
if (-e "$dir/Tree"){

    if (-e "$dir/Tree/3.PhylogeneticTree/rectangular.png"){
        $section->desc("$blank6 本次分析基于<b>核心基因(Core Gene)</b>使用最大似然法构建进化树。");
        # 多图
        my @prefix_tree = qw(rectangular rectangular_bl slanted circular);
        my @fig_tree = map {"./Tree/3.PhylogeneticTree/$_.png"} @prefix_tree;
        my @labels = qw(矩形 矩形真实距离 倾斜 圆形);
        Image($section, \@fig_tree, \@labels, "进化树图");

    }else{
        $section->desc("无进化树结果");
    }

}else{
    $section->desc("未选择进化树选项");
}

################################################################################################################################
#分型
my $section;
$section = $report->section(id => "typing");
$section->menu("分型结果");
if (-e "$dir/Typing"){

    if (-e "$dir/Typing/$sample/mlst.tsv"){
    $section->submenu("$sample mlst分型结果");
    $section->tsv2html(-file=> "$dir/Typing/$sample/mlst.tsv"
            ,-name=> "$sample mlst分型结果"
            ,-header=> 1);

    }else{
        $section->desc("$sample 无mlst分型结果");
    }

    ###############################################################
    if (-e "$dir/Typing/$sample/typing.tsv"){
    $section->submenu("样本 $sample 病原微生物分型");

    $section->tsv2html(-file=> "$dir/Typing/$sample/typing.tsv"
            ,-name=> "$sample 病原微生物分型"
            ,-header=> 1);
    }else{
        $section->desc("$sample 无病原微生物分型结果");
    }

}else{
    $section->desc("未选择分型选项");
}

# 最后
$report->write();







