import logging
from pathlib import Path
from subprocess import run


class ONTFastq:
    """SOALR Nanopore FASTQ 数据处理"""

    def __init__(self, indict, workdir) -> None:
        self.dict_in = indict
        self.dir_work = Path(workdir).resolve()

    def check_empty_files(self, dir_fq: Path) -> list:
        """检查FASTQ目录可能存在的空FASTQ文件, 返回非空FASTQ路径列表"""
        fastqs = []
        for fl in dir_fq.iterdir():
            # 只有一条记录的test.fastq.gz文件大小为171
            if fl.stat().st_size > 200:
                fastqs.append(str(fl))
            assert fastqs, "文件夹下不存在.fastq/.fq/.gz!"
        return fastqs

    def copy_fastq(self):
        logging.info("复制原始数据.")
        self.dir_work.joinpath(".rawdata").mkdir(parents=True, exist_ok=True)
        for smp in self.dict_in['samples']:
            seqdata1 = Path(self.dict_in['samples'][smp]['sequencing_data1'])
            if seqdata1.is_file():
                if seqdata1.suffix != ".gz":
                    cml = f"ln -sf {seqdata1} {self.dir_work}/.rawdata/{smp}.fastq"
                else:
                    cml = f"zcat {seqdata1} > {self.dir_work}/.rawdata/{smp}.fastq"
            elif seqdata1.is_dir():
                fastqs = self.check_empty_files(seqdata1)
                fastqs_sep_by_space = " ".join(fastqs)
                if list(seqdata1.glob('*.fastq')) or list(seqdata1.glob('*.fq')):
                    cml = f"cat {fastqs_sep_by_space} > {self.dir_work}/.rawdata/{smp}.fastq"
                else:
                    cml = f"zcat {fastqs_sep_by_space} > {self.dir_work}/.rawdata/{smp}.fastq"
            else:
                raise TypeError(f"{seqdata1} - 文件类型错误!")
            logging.info("CommandLine: " + cml)
            run(cml, shell=True)


class HYBfastq(ONTFastq):
    """SOALR Nanopore+Illumina FASTQ 数据处理"""

    def __init__(self, indict, workdir) -> None:
        super().__init__(indict, workdir)

    def check_hyb_dir_structure(self, sdata: Path):
        """混合组装目录
        1. 包含 NGS/TGS 目录
        2. NGS 目录下面有2个 FASTQ 文件
        3. TGS 目录下至少有一个非空 FASTQ
        """
        assert {"NGS", "TGS"}.issubset([str(p.name) for p in sdata.iterdir()]), "不存在NGS或TGS目录!"
        # illumina
        ifqs = list(sdata.joinpath("NGS").iterdir())
        assert len(ifqs) == 2, "NGS FASTQ 文件数不对!"
        # nanopore
        nfqs = []
        for fl in sdata.joinpath("TGS").iterdir():
            if fl.stat().st_size > 200:
                nfqs.append(str(fl))
        assert nfqs, "文件夹下不存在.fastq/.fq/.gz!"
        return ifqs, nfqs

    def copy_fastq(self):
        """复制FASTQ. sample.1.fastq, sample.2.fastq, sample.3.fastq"""
        logging.info("复制原始数据.")
        self.dir_work.joinpath(".rawdata").mkdir(parents=True, exist_ok=True)
        for smp in self.dict_in["samples"]:
            seqdata1 = Path(self.dict_in["samples"][smp]["sequencing_data1"]).resolve()
            ilmnfqs, ontfqs = self.check_hyb_dir_structure(seqdata1)
            # illumina
            if ilmnfqs[0].suffix == ".gz":
                cml = f"""
zcat {ilmnfqs[0]} > {self.dir_work}/.rawdata/{smp}.1.fastq
zcat {ilmnfqs[1]} > {self.dir_work}/.rawdata/{smp}.2.fastq
"""
            else:
                cml = f"""
ln -sf {ilmnfqs[0]} {self.dir_work}/.rawdata/{smp}.1.fastq
ln -sf {ilmnfqs[1]} {self.dir_work}/.rawdata/{smp}.2.fastq
"""
            logging.info("CommandLine: " + cml)
            run(cml, shell=True)
            # ont
            fastqs_sep_by_space = " ".join(ontfqs)
            if Path(ontfqs[0]).suffix == ".gz":
                cml = f"zcat {fastqs_sep_by_space} > {self.dir_work}/.rawdata/{smp}.3.fastq"
            else:
                cml = f"cat {fastqs_sep_by_space} > {self.dir_work}/.rawdata/{smp}.3.fastq"
            logging.info("CommandLine: " + cml)
            run(cml, shell=True)

    def generate_samplesheet(self):
        """生成nextflow流程输入samplesheet文件"""
        with open(f"{self.dir_work}/.samplesheet.csv", "w") as g:
            for smp in self.dict_in['samples']:
                g.write(','.join([smp, f'.rawdata/{smp}.1.fastq', f'.rawdata/{smp}.2.fastq',
                                  f'.rawdata/{smp}.3.fastq']) + '\n')
