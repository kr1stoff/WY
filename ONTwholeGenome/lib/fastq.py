import logging
from pathlib import Path
from subprocess import run


class ONTFastq:
    """SOALR Nanopore FASTQ 数据处理"""

    def __init__(self, indict, workdir) -> None:
        self.dict_in = indict
        self.dir_work = workdir

    def check_empty_files(self, dir_fq: Path) -> list:
        """检查FASTQ目录可能存在的空FASTQ文件, 返回非空FASTQ路径列表"""
        fastqs = []
        for fl in dir_fq.iterdir():
            # 只有一条记录的test.fastq.gz文件大小为171
            if fl.stat().st_size > 200:
                fastqs.append(str(fl))
        return fastqs

    def copy_fastq(self):
        logging.info("复制原始数据.")
        for smp in self.dict_in['samples']:
            seqdata1 = Path(self.dict_in['samples'][smp]['sequencing_data1'])
            if seqdata1.is_file():
                if seqdata1.suffix != ".gz":
                    cml = f"ln -sf {seqdata1} {self.dir_work}/.rawdata/{smp}.fastq"
                else:
                    cml = f"zcat {seqdata1} > {self.dir_work}/.rawdata/{smp}.fastq"
            elif seqdata1.is_dir():
                fastqs = self.check_empty_files(seqdata1)
                if not fastqs:
                    raise FileNotFoundError(f"{seqdata1} - 文件夹下不存在.fastq/.fq/.gz!")
                fastqs_sep_by_space = " ".join(fastqs)
                if list(seqdata1.glob('*.fastq')) or list(seqdata1.glob('*.fq')):
                    cml = f"cat {fastqs_sep_by_space} > {self.dir_work}/.rawdata/{smp}.fastq"
                else:
                    cml = f"zcat {fastqs_sep_by_space} > {self.dir_work}/.rawdata/{smp}.fastq"
            elif not seqdata1.exists():
                raise FileNotFoundError(f"{seqdata1} - 文件不存在!")
            else:
                raise TypeError(f"{seqdata1} - 文件类型错误!")
            logging.info("CommandLine: " + cml)
            run(cml, shell=True)
