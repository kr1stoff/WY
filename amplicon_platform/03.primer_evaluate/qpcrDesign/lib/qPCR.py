import logging
import os
import subprocess
from pathlib import Path
import yaml

# 配置文件
d_lib = Path(__file__).absolute().parent
d_base = d_lib.parent
d_config = d_base.joinpath("conf")
f_software = d_config.joinpath("software.yml")
software = yaml.safe_load(open(f_software, 'r').read())
tmpltp3 = d_base.joinpath("etc/template.p3")

# 日志格式
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class Qprimer():
    """
    qPCR引物设计流程.
    """

    def __init__(self, fasta: Path, out: Path, dryrun):
        """
        Init the object

        :param fasta: The fasta sequence
        :param out: The output dir for the result
        :param run: Run the command direct or just generate
        """
        self.d_out = self.set_out(out)
        self.fasta = fasta.absolute()
        self.prepare = self.d_out.joinpath("Prepare")
        self.prepare.mkdir(exist_ok=True)
        self.dryrun = dryrun
        self.bin = d_base.joinpath("bin")

        # 常用软件
        self.python = software["python"]
        self.seqkit = software["seqkit"]

        # 命令
        self.command = []

    def format(self):
        """
        Upcase the fasta file
        """
        f_out = self.prepare.joinpath("sequence.fasta")
        cmd = f"{self.seqkit} seq {self.fasta} -u > {f_out}"
        self.command.append(cmd)
        self.fasta = f_out

    def internal_design(self, size_range='70-110'):
        """
        内引物及探针的设计
        """
        # TODO: 参数优化
        self.d_internal = self.d_out.joinpath("Internal")
        os.makedirs(self.d_internal, exist_ok=True)

        #引物设计
        cmd = f'{self.python} {self.bin}/internal_primer_design.py -f {self.fasta} -o {self.d_internal} -r {size_range}'
        self.command.append(cmd)

        # Primer3 结果整理
        f_out = self.d_internal.joinpath("primers.txt")
        cmd = f"{self.python} {self.bin}/internal_primer_parser.py -f {f_out} --ref {self.fasta} -o {self.d_internal}"
        self.command.append(cmd)

    def set_out(self, out: Path):
        """
        Set the out put dir

        :param out: The out put dir
        """
        d_out = out.absolute()
        d_out.mkdir(exist_ok=True)
        return out

    def run(self):
        """
        Run the generated script
        """
        logging.info("Start to run")
        command = f"bash {self.f_script} >{self.f_log} 2>&1"
        subprocess.run(command, shell=True)

    def finish(self):
        """
        Finish the pipe
        """
        # 生成脚本
        self.f_script = os.path.join(self.d_out, "run.sh")
        self.f_log = os.path.join(self.d_out, "run.log")
        logging.info(f"Write the command to script {self.f_script}")
        with open(self.f_script, 'w') as OUT:
            print("set -eu", file=OUT)
            for command in self.command:
                print(command, file=OUT)

        # 执行脚本
        if not self.dryrun:
            self.run()
