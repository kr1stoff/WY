import os
import sys
import logging
from pathlib import Path
import math

sys.path.append(Path(__file__).resolve())
from .TraceWGS import PhyloWGS


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler)


class PhyloSNP(PhyloWGS):
    """
    溯源进化树 - SNP
    支持 单端/双端/FASTA 三种输入
    [inyaml]    输入YAML文件
    """
    def __init__(self, inyaml):
        super().__init__(inyaml)
        self.assign_software_addition()
        self.assign_baseinfo_addition()
        self.cmds = list()

    def assign_software_addition(self):
        self.bcftools = self.dict_soft['bcftools']
        self.snippy = self.dict_soft['snippy']

    def assign_baseinfo_addition(self):
        self.reference = self.dict_input['reference']
        self.kindom = self.dict_input['kindom']

    def snippy_multi(self):
        """snippy多样本一起跑替代snippy-multi, snippy-multi输出目录不能更改,也没有并行的功能"""
        self.mylog.info("变异检测流程")
        _shell_dir = f"{self.outdir}/Shell/snippy_multi"
        hh = open(f"{self.outdir}/Shell/snippy_multi.sh", "wt", encoding="utf-8", newline="")
        snippy_path = os.path.dirname(self.snippy)
        for name in self.dict_sample:
            gh = open(f"{_shell_dir}/{name}.sh", "wt", encoding="utf-8", newline="")
            if isinstance(self.dict_sample[name], str): # 兼容FASTA样本
                option_input = f"--ctgs {self.dict_sample[name]}"
            elif len(self.dict_sample[name]) == 1: # 样本列表，一个样本SE
                option_input = f"--se {self.dict_sample[name][0]}"
            elif len(self.dict_sample[name]) == 2: # 样本列表，两个样本PE
                option_input = f"--R1 {self.dict_sample[name][0]} --R2 {self.dict_sample[name][1]}"
            else:
                logging.error("样本列表数量错误!")
                raise Exception
            gh.write(f"""
export PATH={snippy_path}:$PATH
{self.snippy} {option_input} --ref {self.reference} --cpus {self.thread} --force \
    --prefix {name} --outdir {self.outdir}/1.Variants/{name}
            """)
            hh.write(f"""
bash {_shell_dir}/{name}.sh > {self.outdir}/logs/1.Variant_{name}.sh.out 2> {self.outdir}/logs/1.Variant_{name}.sh.err
            """)
            gh.close()
        hh.close()
        self.cmds.append(f"""
# 并行跑snippy
{self.python} {self.dir_bin}/kparallel.py -p {self.parallel} {self.outdir}/Shell/snippy_multi.sh \
    > {self.outdir}/logs/snippy_multi.out 2> {self.outdir}/logs/snippy_multi.err
        """)

    def merge_vcf(self):
        self.mylog.info("合并变异文件")
        file_merge_vcf = f"{self.outdir}/2.MergeVCF/vcf_list.txt"
        with open(file_merge_vcf, "wt", encoding="utf-8", newline="") as gh:
            for samp in self.dict_sample.keys():
                gh.write(f"{self.outdir}/1.Variants/{samp}/{samp}.vcf.gz\n")
        self.cmds.append(f"""
# 合并变异文件
{self.bcftools} merge -m snps -f PASS,. --force-samples --output-type v \
    --file-list {file_merge_vcf} -o {self.outdir}/2.MergeVCF/{self.library}.vcf \
    2> {self.outdir}/logs/bcftools_merge.err
        """)

    def build_tree(self):
        self.mylog.info("进化树构建")
        # [220728] vcf2phylip.py 加入--fasta参数,fasttree用.fasta建树
        # [231031] yupeiyu 进化树 -m 参数调整, 样本数 * 0.1
        min_samples_locus = math.ceil(len(self.dict_sample.keys()) * 0.1)
        self.cmds.append(f"""
# 进化树构建
{self.python} {self.dir_bin}/vcf2phylip.py -m {min_samples_locus} --fasta \
    -i {self.outdir}/2.MergeVCF/{self.library}.vcf --output-folder {self.outdir}/3.PhylogeneticTree \
    > {self.outdir}/logs/vcf2phylip.out 2> {self.outdir}/logs/vcf2phylip.err
{self.fasttree} -nt {self.outdir}/3.PhylogeneticTree/{self.library}.min{min_samples_locus}.fasta \
    > {self.outdir}/3.PhylogeneticTree/{self.library}.tre 2> {self.outdir}/logs/fasttree.err
        """)

    def remove_intermedia(self):
        self.mylog.info("删除中间文件")
        self.cmds.append("# 删除中间文件")
        self.cmds.append(f"rm -r {self.outdir}/1.Variants/*/ref*")
    
    def upload(self):
        self.mylog.info("上传数据")
        self.cmds.append(f"""# 上传数据
mkdir {self.outdir}/Upload/source
cp -rf {self.dir_src} {self.outdir}/Upload/source
cp -f {self.outdir}/2.MergeVCF/{self.library}.vcf {self.outdir}/Upload/2.MergeVCF
cp -f {self.outdir}/3.PhylogeneticTree/{self.library}.tre {self.outdir}/Upload/3.PhylogeneticTree
        """)
        for tf in self.tree_figs:
            self.cmds.append(f"""
cp -f {self.outdir}/3.PhylogeneticTree/{tf}.png {self.outdir}/Upload/3.PhylogeneticTree
cp -f {self.outdir}/3.PhylogeneticTree/{tf}.svg {self.outdir}/Upload/3.PhylogeneticTree
            """)

    def execute(self, dryrun:bool):
        """总的执行脚本"""
        logging.info("溯源进化树>全基因组SNP流程: 开始分析!")
        self.make_result_dir(["Shell/snippy_multi", "logs", "1.Variants", "2.MergeVCF", "3.PhylogeneticTree"])
        self.snippy_multi()
        self.merge_vcf()
        self.build_tree()
        self.tree_visualization()
        self.remove_intermedia()
        self.upload()
        self.report('snp')
        self.zip_result()
        self.myrun(dryrun)
        logging.info("溯源进化树>全基因组SNP流程: 分析完成!")
