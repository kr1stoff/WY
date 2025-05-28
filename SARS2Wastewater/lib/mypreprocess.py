import logging
from subprocess import run
from pathlib import Path
from collections import namedtuple
import re
import yaml
import pandas as pd


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class PreProcess():
    def __init__(self, inyaml, analysis, platform) -> None:
        """分析数据前处理."""
        self.inyaml = inyaml
        self.analysis = analysis
        self.platform = platform

    def get_params(self):
        """获取软件/数据库等参数"""
        if self.platform == 'ONT':
            yml = Path(__file__).parents[1].joinpath("config/ont.yaml")
        else:
            yml = Path(__file__).parents[1].joinpath("config/config.yaml")
        self.dict_config = yaml.safe_load(open(yml, "rt"))
        Params = namedtuple("Params", self.dict_config.keys())
        self.params = Params(**self.dict_config)

    def parse_inyaml(self):
        logging.info('解析SOLAR YAML文件.')
        self.dict_inyml = yaml.safe_load(open(self.inyaml))
        self.task_id = self.dict_inyml['task_id']
        self.workdir = f'{self.analysis}/{self.task_id}'

    def make_directories(self):
        logging.info('创建结果目录.')
        cml = f"""
        mkdir -p {self.workdir}/.rawdata {self.workdir}/.material {self.workdir}/.log {self.workdir}/Upload
        mkdir -p {self.workdir}/1.qc {self.workdir}/2.align {self.workdir}/3.muts 
        mkdir -p {self.workdir}/4.consensus {self.workdir}/5.lineage {self.workdir}/6.demix
        """
        logging.info(re.sub(' +', ' ', cml))
        run(cml, shell=True)

    def link_fastq(self):
        logging.info('软链接原始数据.')
        dicsmps = self.dict_inyml['samples']
        for smp in self.dict_inyml['samples']:
            sfx = '.gz' if dicsmps[smp]['sequencing_data1'].endswith('.gz') else ''
            if self.platform == 'ILMN':
                if dicsmps[smp]['sequencing_data2'] != '': #双端
                    cml = f"""
                    ln -sf {dicsmps[smp]['sequencing_data1']} {self.workdir}/.rawdata/{smp}.1.fastq{sfx}            
                    ln -sf {dicsmps[smp]['sequencing_data2']} {self.workdir}/.rawdata/{smp}.2.fastq{sfx}
                    """
                else:
                    cml = f"ln -sf {dicsmps[smp]['sequencing_data1']} {self.workdir}/.rawdata/{smp}.fastq{sfx}"
            else: #ONT
                if Path(dicsmps[smp]['sequencing_data1']).is_file(): #单个fastq
                    cml = f"ln -sf {dicsmps[smp]['sequencing_data1']} {self.workdir}/.rawdata/{smp}.fastq{sfx}"
                elif Path(dicsmps[smp]['sequencing_data1']).is_dir():
                    if list(Path(dicsmps[smp]['sequencing_data1']).glob('*.fastq')): #存在.fastq
                        cml = f"cat {dicsmps[smp]['sequencing_data1']}/*.fastq > {self.workdir}/.rawdata/{smp}.fastq"
                    elif list(Path(dicsmps[smp]['sequencing_data1']).glob('*.fq')): #存在.fq
                        cml = f"cat {dicsmps[smp]['sequencing_data1']}/*.fq > {self.workdir}/.rawdata/{smp}.fastq"
                    elif list(Path(dicsmps[smp]['sequencing_data1']).glob('*.fastq.gz')): #存在.fastq
                        cml = f"zcat {dicsmps[smp]['sequencing_data1']}/*.fastq.gz > {self.workdir}/.rawdata/{smp}.fastq"
                    elif list(Path(dicsmps[smp]['sequencing_data1']).glob('*.fq.gz')): #存在.fq
                        cml = f"zcat {dicsmps[smp]['sequencing_data1']}/*.fq.gz > {self.workdir}/.rawdata/{smp}.fastq"
                    else:
                        raise Exception('输入FASTQ文件夹为空, 或不存在FASTQ文件!')
            logging.info(re.sub(' +', ' ', cml))
            run(cml, shell=True)

    def generate_snakemake_config(self):
        logging.info('生成Snakemake配置文件.')
        dicsmkcfg = {}
        dicsmkcfg['workdir'] = str(Path(self.workdir).resolve())
        dicsmkcfg['samples'] = list(self.dict_inyml['samples'].keys())
        dicsmkcfg.update(self.dict_config)
        #引物
        prmseq = self.dict_inyml['primer_seq']
        prmbed = f"{dicsmkcfg['workdir']}/primer.bed"
        #1.有提交文件; 2.文件存在; 3.文件后缀.xlsx
        if prmseq and Path(prmseq).resolve().is_file() and prmseq.endswith(".xlsx"):
            df = pd.read_excel(prmseq, header=2)
            df.to_csv(prmbed, sep='\t', index=False)
            dicsmkcfg['database']['primer_bed'] = prmbed
        with open(f'{self.workdir}/.material/configfile.yaml', 'wt') as g:
            g.write(yaml.safe_dump(dicsmkcfg))
