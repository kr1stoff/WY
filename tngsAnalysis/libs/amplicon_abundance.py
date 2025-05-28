import logging
from pathlib import Path
from subprocess import run
from modules.config import config


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class AmpliconAbundance():
    def __init__(self, fastq1, fastq2, reference, bed, outdir, prefix, dryrun=False):
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.reference = reference
        self.bed = bed
        self.outdir = outdir
        self.prefix = prefix
        self.dryrun = dryrun
        self.myconfig()
        self.cmds = []

    def myconfig(self):
        self.params = config()
        self.threads = self.params.threads
        self.dir_conf = Path(__file__).parents[1].joinpath('conf')
        self.dir_bin = Path(__file__).parents[1].joinpath('bin')

    def mymkdir(self):
        self.cmds.append(f"""
#创建目录
mkdir -p {self.outdir}/tmp
mkdir -p {self.outdir}/1.qc {self.outdir}/1.qc/before {self.outdir}/1.qc/after
mkdir -p {self.outdir}/2.map
        """)

    def link_data(self):
        """软连接原始文件, 相当于格式化文件名了"""
        suffix = 'fq.gz' if Path(self.fastq1).suffix == '.gz' else 'fq'
        self.link_fq1 = f'{self.outdir}/tmp/{self.prefix}.1.{suffix}'
        self.link_fq2 = f'{self.outdir}/tmp/{self.prefix}.2.{suffix}'
        cml_link_fastq2 = ''
        if self.fastq2:
            cml_link_fastq2 = f'ln -sf {Path(self.fastq2).resolve()} {self.link_fq2}'
        self.cmds.append(f"""
#软连接原始数据
ln -sf {Path(self.fastq1).resolve()} {self.link_fq1}
{cml_link_fastq2}
        """)

    def qc(self):
        """质控过滤"""
        if self.fastq2:
            fastq2 = self.link_fq2
            fastp2I = f'-I {self.link_fq2}'
            fastp2O = f'-O {self.outdir}/1.qc/{self.prefix}.clean.2.fq'
            parse_fastp_mode = 'PE'
        else:
            fastq2,fastp2I,fastp2O = '', '', ''
            parse_fastp_mode = 'SE'
        self.cmds.append(f"""
#原始数据QC
{self.params.fastqc} -t {self.threads} --extract -o {self.outdir}/1.qc/before {self.link_fq1} {fastq2} &
{self.params.fastp} --thread {self.threads} {self.params.options_fastp} \
    --adapter_fasta {self.dir_conf}/adapter.fa \
    -j {self.outdir}/1.qc/{self.prefix}.json -h {self.outdir}/1.qc/{self.prefix}.html \
    -i {self.link_fq1} {fastp2I} \
    -o {self.outdir}/1.qc/{self.prefix}.clean.1.fq {fastp2O}
python3 {self.dir_bin}/parse_fastp.json.py \
    -m {parse_fastp_mode} -id {self.prefix} \
    -i {self.outdir}/1.qc/{self.prefix}.json
{self.params.fastqc} -t {self.threads} --extract \
    {self.outdir}/1.qc/{self.prefix}.clean.1.fq {fastp2O.replace('-O ','')} \
    -o {self.outdir}/1.qc/after 
        """)

    def map(self):
        """比对"""
        self.ref = f'{self.outdir}/tmp/{Path(self.reference).name}'
        fastq2 = f'{self.outdir}/1.qc/{self.prefix}.clean.2.fq' if self.fastq2 else ''
        self.cmds.append(f"""
#比对
ln -sf {Path(self.reference).resolve()} {self.ref}
{self.params.bwa} index {self.ref}
{self.params.bwa} mem -t {self.threads} -M -Y -R '@RG\\tID:{self.prefix}\\tSM:{self.prefix}' \
    {self.ref} {self.outdir}/1.qc/{self.prefix}.clean.1.fq {fastq2} \
    > {self.outdir}/2.map/{self.prefix}.raw.sam
#统计过滤
{self.params.python} {self.dir_bin}/count_reads.py --identity 0.8 \
    {self.outdir}/2.map/{self.prefix}.raw.sam \
    -o {self.outdir}/2.map/{self.prefix}.filtered.sam \
    -s {self.outdir}/2.map/{self.prefix}.count_reads.stat
{self.params.samtools} sort -@ {self.threads} {self.outdir}/2.map/{self.prefix}.filtered.sam -O BAM \
    > {self.outdir}/2.map/{self.prefix}.filtered.sort.bam
ln -sf {self.prefix}.filtered.sort.bam {self.outdir}/2.map/{self.prefix}.bam
{self.params.samtools} index {self.outdir}/2.map/{self.prefix}.bam
        """)

    def abundance(self):
        """扩增子丰度计算"""

    def write_and_run(self):
        """命令行写入脚本并"""
        Path(self.outdir).mkdir(exist_ok=True, parents=True)
        with open(f'{self.outdir}/pipe.sh', 'wt', encoding='utf-8', newline='') as g:
            for cml in self.cmds:
                g.write(cml + '\n')
        cmd = f'bash {self.outdir}/pipe.sh &> {self.outdir}/pipe.sh.log'
        if not self.dryrun:
            run(cmd, shell=True)

    def execute(self):
        self.mymkdir()
        self.link_data()
        self.qc()
        self.map()
        self.write_and_run()
