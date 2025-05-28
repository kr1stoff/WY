import sys
from pathlib import Path
import logging
import yaml
import re
from subprocess import run
sys.path.append(Path(__file__).resolve())
from . import general


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler)


class PhyloWGS():
    """
    溯源进化树 - WGS
    [inyaml]    输入YAML文件
    """
    def __init__(self, inyaml):
        self.inyaml = inyaml
        self.locate_home()
        self.assign_software()
        self.assign_parameters()
        self.assign_baseinfo()
        self.mylog = general.MyLoggingInfo()
        self.cmds = list() # 总脚本中的所有命令列表
        
    def locate_home(self):
        """定位文件目录"""
        self.home = Path(__file__).resolve().parents[1]
        self.dir_bin = self.home.joinpath("bin")
        self.dir_conf = self.home.joinpath("conf")
        self.dir_src = self.home.joinpath("src")

    def assign_software(self):
        """配置软件路径"""
        logging.info("配置软件路径")
        self.dict_soft = yaml.safe_load(open(f"{self.dir_conf}/software.yml"))
        self.python = self.dict_soft["python"]
        self.fasttree = self.dict_soft["fasttree"]
        self.Rscript = self.dict_soft["Rscript"]
        self.magick = self.dict_soft["magick"]
        self.perl = self.dict_soft["perl"]
        self.mafft = self.dict_soft['mafft']

    def assign_parameters(self):
        """软件配置相关类参数"""
        logging.info("读取软件参数配置文件")
        self.tree_figs = ["rectangular", "rectangular_bl", "slanted", "circular"]
        self.dict_params = yaml.safe_load(open(f"{self.dir_conf}/parameters.yml"))
        self.thread = self.dict_params["thread"]
        self.parallel = self.dict_params["parallel"]

    def assign_baseinfo(self):
        """基础信息"""
        logging.info("读取输入配置文件")
        self.dict_input = yaml.safe_load(open(self.inyaml))
        self.library = self.dict_input["library"]
        self.outdir = Path(self.dict_input["result_dir"]).joinpath(self.library)
        self.dict_sample = self.dict_input['samples']

    def make_result_dir(self, dirs:list):
        """创建结果目录,Upload在创建一遍"""
        self.mylog.info("创建结果目录")
        # 数字开头的目录才上传到Upload目录
        uploads_release = [dd for dd in dirs if re.match("^\d\.", dd)]
        uploads = list(map(lambda x:f"Upload/{x}", uploads_release))
        general.mymakedirs(dirs+uploads, self.outdir)

    def merge_fasta(self):
        self.mylog.info("合并FASTA文件")
        self.fasta_merged = f"{self.outdir}/1.MergeFA/{self.library}.fa"
        general.merge_wgs_fasta(mfasta=self.fasta_merged, dict_sample=self.dict_sample)

    def msa(self, infasta):
        self.mylog.info("MSA 多序列比对")
        self.cmds.append(
            f"""
# MSA 多序列比对
{self.mafft} --thread {self.thread} --auto --maxiterate 1000 {infasta} \
    > {self.outdir}/2.MSA/{self.library}.aln.fa 2> {self.outdir}/logs/mafft.err
            """)

    def build_tree(self):
        self.mylog.info("进化树构建")
        self.cmds.append(f"""
# 进化树构建
{self.fasttree} -nt {self.outdir}/2.MSA/{self.library}.aln.fa \
    > {self.outdir}/3.PhylogeneticTree/{self.library}.tre 2> {self.outdir}/logs/fasttree.err
        """)

    def tree_visualization(self):
        self.mylog.info("进化树可视化")
        self.cmds.append(f"""
# 进化树可视化
{self.Rscript} {self.dir_bin}/tree.R {self.outdir}/3.PhylogeneticTree/{self.library}.tre \
    > {self.outdir}/logs/rtree.out 2> {self.outdir}/logs/rtree.err
{self.python} {self.dir_bin}/magick.py {self.magick} {self.outdir}/3.PhylogeneticTree \
    > {self.outdir}/logs/magick.out 2> {self.outdir}/logs/magick.err
        """)
    
    def upload(self):
        self.mylog.info("上传数据")
        cml = f"""# 上传数据
mkdir {self.outdir}/Upload/source
cp -rf {self.dir_src} {self.outdir}/Upload/source
cp -f {self.fasta_merged} {self.outdir}/Upload/1.MergeFA
cp -f {self.outdir}/2.MSA/{self.library}.aln.fa {self.outdir}/Upload/2.MSA
cp -f {self.outdir}/3.PhylogeneticTree/{self.library}.tre {self.outdir}/Upload/3.PhylogeneticTree
        """
        for tf in self.tree_figs:
            cml += f"""
cp -f {self.outdir}/3.PhylogeneticTree/{tf}.png {self.outdir}/Upload/3.PhylogeneticTree
cp -f {self.outdir}/3.PhylogeneticTree/{tf}.svg {self.outdir}/Upload/3.PhylogeneticTree
        """
        self.cmds.append(cml)

    def report(self, pipe='wgs'):
        """pipe  ::  {wgs, snp, core}"""
        self.mylog.info("生成报告")
        self.cmds.append(f"""
# 生成报告
{self.perl} {self.dir_bin}/report.pl {self.outdir}/Upload {pipe} \
    > {self.outdir}/logs/report.out 2> {self.outdir}/logs/report.err
                        """)

    def zip_result(self):
        self.mylog.info("结果压缩包")
        self.cmds.append(f"""
# 结果压缩包
cd {self.outdir}
zip -r {self.library}.zip Upload
        """)

    def myrun(self, dryrun:bool):
        """运行步骤"""
        self.mylog.info('运行总脚本')
        all_script = f'{self.outdir}/Shell/all.sh'
        with open(all_script, 'wt', encoding='utf-8', newline='') as gh:
            for cmd in self.cmds:
                gh.write(cmd+'\n')
        if not dryrun:
            run(f'bash {all_script} > {self.outdir}/logs/all.sh.out 2> {self.outdir}/logs/all.sh.err', shell=True)

    def execute(self, dryrun:bool):
        logging.info('溯源进化树>全基因组多序列比对流程: 开始分析!')
        self.make_result_dir(['1.MergeFA', '2.MSA', '3.PhylogeneticTree', 'logs', 'Shell'])
        self.merge_fasta()
        self.msa(self.fasta_merged)
        self.build_tree()
        self.tree_visualization()
        self.upload()
        self.report()
        self.zip_result()
        self.myrun(dryrun)
        logging.info('溯源进化树>全基因组多序列比对流程: 分析完成!')
