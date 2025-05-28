import logging
from subprocess import run
from pathlib import Path
import os
import re


class Analysis():
    def __init__(self, prprcs) -> None:
        """分析"""
        self.params = prprcs.params
        self.workdir = prprcs.workdir
        self.dict_inyml = prprcs.dict_inyml

    def run_snakemake_pipe(self) -> None:
        """运行snakemake流程. prprcs是前处理的信息对象."""
        logging.info('运行snakemake流程.')
        smkenv = Path(self.params.path['snakemake']).resolve().parents[1]
        activate = Path(self.params.path['activate']).resolve()
        smkfl = Path(self.params.custom['wastewater_smkpipe']).resolve()
        num_cpu = os.cpu_count()
        configfile = f'{self.workdir}/.material/configfile.yaml'
        cml = f"""
        source {activate} {smkenv}
        snakemake -c {num_cpu} -s {smkfl} --configfile {configfile}
        """
        logging.info(re.sub(' +', ' ', cml))
        run(cml, shell=True, executable='/bin/bash', capture_output=True)  # 接收输出防止阻塞(NewBing)

    def handle_user_predict_lineages(self):
        """
        [230418]用户预测分型
        SOLAR前端-更多样本信息-样本备注
        规则: 逗号分隔分型, 例 AY.100,XBB.2.6,BF.7.15,Q.1
        MySQL中字段名称为'sample_properties'
        """
        logging.info('处理用户预测分型.')
        for smp in self.dict_inyml['samples']:
            usrlngstr = self.dict_inyml['samples'][smp]['sample_properties']
            allabdctbl = f'{self.workdir}/6.demix/{smp}.abundance.txt'
            usrabdctbl = f'{self.workdir}/6.demix/{smp}.abundance_user_predict.txt'
            if usrlngstr:  # 如果填了备注就生成用户预测分型表
                usrlngstr = usrlngstr.replace('，', ',')  # 预防输入中文逗号
                usrlngs = [l for l in usrlngstr.strip().split(',') if l]
                with open(allabdctbl) as f, open(usrabdctbl, 'w') as g:
                    g.write(next(f))
                    for line in f:
                        if line.split('\t')[0] in usrlngs:
                            g.write(line)
