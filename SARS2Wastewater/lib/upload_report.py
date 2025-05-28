import logging
from subprocess import run
from pathlib import Path
import re


class UploadReport():
    def __init__(self, prprcs) -> None:
        """上传&生成报告"""
        self.params = prprcs.params
        self.workdir = prprcs.workdir
        self.dict_inyml = prprcs.dict_inyml
        self.platform = prprcs.platform

    def upload(self):
        logging.info('复制Upload结果文件.')
        dricn = Path(__file__).resolve().parents[1].joinpath('icons')
        for smp in self.dict_inyml['samples']:
            cml = f"""
            cd {self.workdir}
            #quality control
            mkdir -p Upload/{smp}/1.qc
            cp -rf 1.qc/{smp}/before 1.qc/{smp}/after \
                1.qc/{smp}.html 1.qc/{smp}.basic.stat.txt \
                -t Upload/{smp}/1.qc
            #alignment
            mkdir -p Upload/{smp}/2.align
            cp -rf 2.align/{smp}.trmsrt.bam* Upload/{smp}/2.align
            #demix & coverage
            mkdir -p Upload/{smp}/6.demix
            cp -rf 6.demix/{smp}.coverage.png 6.demix/{smp}.coverage.txt \
                6.demix/{smp}.abundance.txt 6.demix/{smp}.pie.png \
                -t Upload/{smp}/6.demix
            #mutations
            mkdir -p Upload/{smp}/3.muts
            cp -rf 3.muts/{smp}.display.tsv Upload/{smp}/3.muts
            #consensus
            mkdir -p Upload/{smp}/4.consensus
            cp -rf 4.consensus/{smp}.consensus.fa Upload/{smp}/4.consensus
            #icons
            mkdir Upload/{smp}/source
            cp -rf {dricn} Upload/{smp}/source
            """
            #双端|单端
            if self.dict_inyml['samples'][smp]['sequencing_data2'] != '': 
                cml += f"""
                cp -rf 1.qc/{smp}/before/{smp}.1_fastqc/Images/per_base_quality.png Upload/{smp}/1.qc/{smp}.1_before_per_base_quality.png
                cp -rf 1.qc/{smp}/after/{smp}.clean.1_fastqc/Images/per_base_quality.png Upload/{smp}/1.qc/{smp}.1_after_per_base_quality.png
                cp -rf 1.qc/{smp}/before/{smp}.2_fastqc/Images/per_base_quality.png Upload/{smp}/1.qc/{smp}.2_before_per_base_quality.png
                cp -rf 1.qc/{smp}/after/{smp}.clean.2_fastqc/Images/per_base_quality.png Upload/{smp}/1.qc/{smp}.2_after_per_base_quality.png
                """
            else:
                cml += f"""
                cp -rf 1.qc/{smp}/before/{smp}_fastqc/Images/per_base_quality.png Upload/{smp}/1.qc/{smp}_before_per_base_quality.png
                cp -rf 1.qc/{smp}/after/{smp}.clean_fastqc/Images/per_base_quality.png Upload/{smp}/1.qc/{smp}_after_per_base_quality.png
                """
            #用户预测分型表
            if Path(f'{self.workdir}/6.demix/{smp}.abundance_user_predict.txt').exists():
                cml += f'cp -rf 6.demix/{smp}.abundance_user_predict.txt -t Upload/{smp}/6.demix'
            logging.info(re.sub(' +', ' ', cml))
            run(cml, shell=True, executable='/bin/bash', capture_output=True)

    def uploadONT(self):
        logging.info('复制UploadONT结果文件.')
        dricn = Path(__file__).resolve().parents[1].joinpath('icons')
        for smp in self.dict_inyml['samples']:
            cml = f"""
            cd {self.workdir}
            #quality control
            mkdir -p Upload/{smp}/1.qc
            cp -rf 1.qc/{smp}.cln.stats.png 1.qc/{smp}.raw.stats.png 1.qc/{smp}.stat.txt -t Upload/{smp}/1.qc
            #alignment
            mkdir -p Upload/{smp}/2.align
            cp -rf 2.align/{smp}.trmst.bam* -t Upload/{smp}/2.align
            #demix & coverage
            mkdir -p Upload/{smp}/6.demix
            cp -rf 6.demix/{smp}.abundance.txt 6.demix/{smp}.coverage.png 6.demix/{smp}.coverage.txt 6.demix/{smp}.pie.png \
                -t Upload/{smp}/6.demix
            #mutations
            mkdir -p Upload/{smp}/3.muts
            cp -rf 3.muts/{smp}.display.tsv -t Upload/{smp}/3.muts
            #consensus
            mkdir -p Upload/{smp}/4.consensus
            cp -rf 4.consensus/{smp}.consensus.fa Upload/{smp}/4.consensus
            #icons
            mkdir Upload/{smp}/source
            cp -rf {dricn} Upload/{smp}/source
            """
            #用户预测分型表
            if Path(f'{self.workdir}/6.demix/{smp}.abundance_user_predict.txt').exists():
                cml += f'cp -rf 6.demix/{smp}.abundance_user_predict.txt -t Upload/{smp}/6.demix'
            logging.info(re.sub(' +', ' ', cml))
            run(cml, shell=True, executable='/bin/bash', capture_output=True)

    def report(self):
        logging.info('生成HTML报告.')
        if self.platform == 'ILMN':
            plrpt = Path(__file__).resolve().parents[1].joinpath('bin/report.pl')
        else:
            plrpt = Path(__file__).resolve().parents[1].joinpath('bin/reportONT.pl')
        for smp in self.dict_inyml['samples']:
            cml = f'{self.params.path["perl"]} {plrpt} {self.workdir}/Upload/{smp} {smp}'
            logging.info(cml)
            run(cml, shell=True, executable='/bin/bash', capture_output=True)

    def zip_upload(self):
        logging.info('压缩Upload结果文件夹.')
        for smp in self.dict_inyml['samples']:
            cml = f"""
            cd {self.workdir}/Upload
            zip -r {smp}.zip {smp}
            """
            logging.info(re.sub(' +', ' ', cml))
            run(cml, shell=True, executable='/bin/bash', capture_output=True)
