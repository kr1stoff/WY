import os
import sys
import logging
from pathlib import Path
sys.path.append(Path(__file__).resolve())
from .TraceWGS import PhyloWGS
from . import general
import pdb


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler)


class PhyloCORE(PhyloWGS):
    """
    溯源进化树 - CORE
    [inyaml]    输入YAML文件
    """
    def __init__(self, inyaml):
        super().__init__(inyaml)
        self.assign_baseinfo_addition()
        self.assign_software_addition()

    def assign_baseinfo_addition(self):
        self.kindom = self.dict_input['kindom']
    
    def assign_software_addition(self):
        self.prodigal = self.dict_soft['prodigal']
        self.genemark = self.dict_soft['gmhmmp']
        self.seqtk = self.dict_soft['seqtk']
        self.orthofinder = self.dict_soft['orthofinder']

    def format_fasta_multi(self):
        self.mylog.info('格式化输入')
        self.dict_sample_formatted = dict()
        for name, path in self.dict_sample.items():
            outfile = f'{self.outdir}/FormatFASTA/{name}.fna'
            general.format_fasta(name, path, outfile)
            self.dict_sample_formatted[name] = outfile

    def run_prodigal(self):
        """gene_predict: 细菌病毒用prodigal并行"""
        self.mylog.info('prodigal基因预测')
        _shell_dir = f'{self.outdir}/Shell/prodigal_multi'
        os.makedirs(_shell_dir, exist_ok=True)
        hh = open(f'{self.outdir}/Shell/prodigal_multi.sh', 'wt', encoding='utf-8', newline='')
        for name, fasta in self.dict_sample_formatted.items():
            gh = open(f'{_shell_dir}/{name}.sh', 'wt', encoding='utf-8', newline='')
            gh.write(f"""
# export PATH={Path(self.prodigal).parent}:$PATH
{self.prodigal} -i {fasta} \
    -a {self.outdir}/1.Orthologue/GenePrediction/{name}.faa \
    -o {self.outdir}/1.Orthologue/GenePrediction/{name}.genes
            """)
            hh.write(f"""
bash {_shell_dir}/{name}.sh \
    > {self.outdir}/logs/prodigal_{name}.out 2> {self.outdir}/logs/prodigal_{name}.err
            """)
            gh.close()
        hh.close()
        self.cmds.append(f"""
# # prodigal 基因预测
{self.python} {self.dir_bin}/kparallel.py -p {self.parallel} {self.outdir}/Shell/prodigal_multi.sh \
    > {self.outdir}/logs/prodigal_multi.out 2> {self.outdir}/logs/prodigal_multi.err
        """)

    def run_genemark(self):
        """gene_predict: 真菌用genemark并行"""
        self.mylog.info('genemark基因预测')
        awk_options = """'{if ($0~/^>/) {split($0,a,\"\\t>\");printf \">%s_%s\\n\",a[2],a[1]} else {print $0}}'"""
        _shell_dir = f'{self.outdir}/Shell/genemark_multi'
        os.makedirs(_shell_dir, exist_ok=True)
        hh = open(f'{self.outdir}/Shell/genemark_multi.sh', 'wt', encoding='utf-8', newline='')
        for name, fasta in self.dict_sample_formatted.items():
            gh = open(f'{_shell_dir}/{name}.sh', 'wt', encoding='utf-8', newline='')
            gh.write(f"""
# export PATH={Path(self.genemark).parent}:$PATH
{self.genemark} -a -d -f G -m {self.dict_params['MetaGeneMarkMod']} {fasta} \
    -A {self.outdir}/1.Orthologue/GenePrediction/{name}.raw.faa \
    -o {self.outdir}/1.Orthologue/GenePrediction/{name}.genes \
    -D {self.outdir}/1.Orthologue/GenePrediction/{name}.fna
awk {awk_options} {self.outdir}/1.Orthologue/GenePrediction/{name}.raw.faa \
    > {self.outdir}/1.Orthologue/GenePrediction/{name}.faa
            """)
            hh.write(f"""
bash {_shell_dir}/{name}.sh > {self.outdir}/logs/genemark_{name}.out 2> {self.outdir}/logs/genemark_{name}.err
            """)
            gh.close()
        hh.close()
        self.cmds.append(f"""
# genemark 基因预测
{self.python} {self.dir_bin}/kparallel.py -p {self.parallel} {self.outdir}/Shell/genemark_multi.sh \
    > {self.outdir}/logs/genemark_multi.out 2> {self.outdir}/logs/genemark_multi.err
        """)
        # pdb.set_trace()

    def genepred_lenfilter(self):
        """预测基因长度筛选，蛋白长度>33bp(核酸长度>100bp)"""
        self.mylog.info('预测基因长度筛选')
        _shell_dir = f'{self.outdir}/Shell/seqtkseq_multi'
        os.makedirs(_shell_dir, exist_ok=True)
        hh = open(f'{self.outdir}/Shell/seqtkseq_multi.sh', 'wt', encoding='utf-8', newline='')
        for name in self.dict_sample:
            gh = open(f'{_shell_dir}/{name}.sh', 'wt', encoding='utf-8', newline='')
            gh.write(f"""
{self.seqtk} seq -L 33 {self.outdir}/1.Orthologue/GenePrediction/{name}.faa \
    > {self.outdir}/1.Orthologue/GenePrediction/{name}_L33.faa
cp -rf {self.outdir}/1.Orthologue/GenePrediction/{name}_L33.faa {self.outdir}/1.Orthologue/orthofinder_indir/{name}_L33.faa
            """)
            hh.write(f"""
bash {_shell_dir}/{name}.sh > {self.outdir}/logs/seqtkseq_{name}.out 2> {self.outdir}/logs/seqtkseq_{name}.err
            """)
            gh.close()
        hh.close()
        self.cmds.append(f"""
# 预测基因筛选
{self.python} {self.dir_bin}/kparallel.py -p {self.parallel} {self.outdir}/Shell/seqtkseq_multi.sh \
    > {self.outdir}/logs/seqtkseq_multi.out 2> {self.outdir}/logs/seqtkseq_multi.err
        """)

    def gene_predict(self):
        """分病毒细菌和真菌两种方向"""
        self.mylog.info("基因预测")
        os.makedirs(f"{self.outdir}/1.Orthologue/GenePrediction", exist_ok=True)
        os.makedirs(f"{self.outdir}/1.Orthologue/orthofinder_indir", exist_ok=True)
        if self.kindom in ["bacteria", "viruses"]:
            self.run_prodigal()
        elif self.kindom == "fungi":
            self.run_genemark()
        else:
            logging.error("只允许 <bacteria/viruses/fungi>")
            raise Exception("错误终止分析!")
        # 长度筛选
        self.genepred_lenfilter()        

    def run_orthofinder(self):
        """orthologue同源基因"""
        with open(f'{self.outdir}/Shell/run_orthofinder.sh', 'wt', encoding='utf-8', newline='') as gh:
            gh.write(f"""
ulimit -n 99999 #调高打开文件数限制
export PATH={Path(self.orthofinder).parent}:$PATH
{self.orthofinder} -t {self.thread} -S diamond -f {self.outdir}/1.Orthologue/orthofinder_indir \
    -n {self.library} -o {self.outdir}/1.Orthologue/orthofinder_out
        """)
        self.cmds.append(f"""
# 同源基因
bash {self.outdir}/Shell/run_orthofinder.sh > {self.outdir}/logs/run_orthofinder.out 2> {self.outdir}/logs/run_orthofinder.err
        """)

    def orthologue(self):
        """1.找同源基因 2.单拷贝基因多序列比对 3.合并成长基因序列"""
        self.mylog.info("核心基因")
        self.file_supergene = f"{self.outdir}/1.Orthologue/{self.library}_supergene.faa"
        # 同源基因
        self.run_orthofinder()
        self.cmds.append(f"""
# 单拷贝基因集
{self.python} {self.dir_bin}/single_copy_orthologue_sequences.py \
    {self.outdir}/1.Orthologue/orthofinder_out/Results_{self.library}/Single_Copy_Orthologue_Sequences \
    {self.inyaml} {self.file_supergene} \
    > {self.outdir}/logs/single_copy_orthologue_sequences.out 2> {self.outdir}/logs/single_copy_orthologue_sequences.err
        """)

    def upload(self):
        self.mylog.info("上传数据")
        cml = f"""# 上传数据
mkdir {self.outdir}/Upload/source
cp -rf {self.dir_src} {self.outdir}/Upload/source
cp -f {self.file_supergene} {self.outdir}/Upload/1.Orthologue
cp -f {self.outdir}/2.MSA/{self.library}.aln.fa {self.outdir}/Upload/2.MSA
cp -f {self.outdir}/3.PhylogeneticTree/{self.library}.tre {self.outdir}/Upload/3.PhylogeneticTree
        """
        for tf in self.tree_figs:
            cml += f"""
cp -f {self.outdir}/3.PhylogeneticTree/{tf}.png {self.outdir}/Upload/3.PhylogeneticTree
cp -f {self.outdir}/3.PhylogeneticTree/{tf}.svg {self.outdir}/Upload/3.PhylogeneticTree
            """
        self.cmds.append(cml)

    def execute(self, dryrun:bool):
        logging.info("溯源进化树>核心基因流程: 开始分析!")
        self.make_result_dir(["1.Orthologue", "2.MSA", "3.PhylogeneticTree", "FormatFASTA", "logs", "Shell"])
        self.format_fasta_multi()
        self.gene_predict()
        self.orthologue()
        self.msa(self.file_supergene)
        self.build_tree()
        self.tree_visualization()
        self.upload()
        self.report('core')
        self.zip_result()
        self.myrun(dryrun)
        logging.info("溯源进化树>核心基因流程: 分析完成!")
