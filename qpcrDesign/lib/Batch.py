import re
import logging
from subprocess import run
from pathlib import Path
from collections import namedtuple
import yaml

# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class BatchQPCR():
    def __init__(self, mypath, reference, bed, out) -> None:
        self.mypath = mypath
        self.ref = Path(reference).resolve()
        self.bed = Path(bed).resolve()
        self.out = Path(out).resolve()
        self.get_params()

    def get_params(self):
        """获取软件路径"""
        yml = Path(__file__).parents[1].joinpath("conf/software.yml")
        self.dict_conf = yaml.safe_load(open(yml, "rt"))
        Params = namedtuple("Params", self.dict_conf.keys())
        self.params = Params(**self.dict_conf)

    def get_region_fasta(self):
        """从BED获取FASTA, 并拆分."""
        cmd = f"""
        set -eu
        cd {self.out}
        {self.params.bedtools} getfasta -fi {self.ref} -bed {self.bed} -fo region.fa
        {self.params.seqkit} split2 -p $(grep '>' -c region.fa) region.fa
        """
        self.myrun(cmd)

    def design1(self):
        """primer3引物设计, 初次设计使用70-110扩增产物长度阈值."""
        self.regprt = [] #region.part_*
        with open(self.out.joinpath('batch70.sh'), 'w') as g:
            for fa in self.out.glob('region.fa.split/region.part_*.fa'):
                cmd = f"""
                {self.params.python} {self.mypath} qpcr -o {self.out}/{fa.stem} -f {fa}
                """
                g.write(self.format_command(cmd) + '\n')
                self.regprt.append(str(fa.stem))
        cmd = f"cat {self.out.joinpath('batch70.sh')} | {self.params.parallel} -j 8"
        self.myrun(cmd)

    def design2(self):
        """primer3引物设计, 针对70-110阈值不成功的区域再次设计,使用50-150扩增产物长度阈值."""
        with open(self.out.joinpath('batch50.sh'), 'w') as g:
            for rp in self.regprt:
                if len(open(self.out.joinpath(f'{rp}/Internal/primer.xls')).readlines()) > 1: #有结果
                    continue
                cmd = f"""
                {self.params.python} {self.mypath} qpcr -r 50-150 -o {self.out}/{rp} -f {self.out}/region.fa.split/{rp}.fa
                """
                g.write(self.format_command(cmd) + '\n')
            cmd = f"cat {self.out.joinpath('batch50.sh')} | {self.params.parallel} -j 8"
            self.myrun(cmd)

    def get_exp_txt(self):
        """合并结果, 根据BED添加NAME列."""
        #读BED,获取region和name字典
        dicregnm = {}
        with open(self.bed) as f:
            for line in f:
                llst = line.strip().split('\t')
                #{(NZ_CP015941.1,978060,978813): Legionella pneumophila-group_12296, ...}
                dicregnm[(llst[0], llst[1], llst[2])] = llst[3] #和'BED NAME'完全一致
        #写exp.txt
        with open(self.out.joinpath('exp.tmp.txt'), 'w') as g:
            for rp in self.regprt:
                num = 1 #局部计数
                with open(self.out.joinpath(f'{rp}/Internal/primer.xls')) as f:
                    lines = f.readlines()
                if len(lines) == 1:
                    continue #没结果
                head = lines[0]
                head = head.strip() + '\tNAME\n'
                g.write(head)
                for line in lines[1:]:
                    llst = line.strip().split('\t')
                    loc3 = tuple(llst[-3:]) #chrom,start,end
                    name = dicregnm[loc3]
                    g.write(line.strip() + f'\t{name}_{num}\n')
                    num += 1
        #230901去重exp.txt. col 7,13,19 左引物,右引物,探针
        cmd = f"{self.params.csvtk} uniq -tH -f 7,13,19 {self.out.joinpath('exp.tmp.txt')} > {self.out.joinpath('exp.txt')}"
        self.myrun(cmd)

    def myrun(self, cmd):
        """run固定参数"""
        logging.info(self.format_command(cmd))
        res = run(cmd, shell=True, encoding='utf-8', executable='/bin/bash', capture_output=True)
        with open(f'{self.out}/batch.log', 'a') as g:
            g.write(res.stdout + res.stderr)

    def format_command(self, cmd):
        """格式化命令"""
        fmtcmd = re.sub(' +', ' ', cmd.strip())
        return fmtcmd
