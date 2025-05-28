#!/usr/bin/env python

from pathlib import Path
import logging
import yaml
from subprocess import run
import re
from multiprocessing import cpu_count


class Pipe:
    """ONT+ILMN混装细菌全基因组分析流程"""

    def __init__(self, inyaml, analysis, conda_activate) -> None:
        self.yaml_in = inyaml
        self.dir_proj = Path(__file__).resolve().parents[1]
        self.conda_activate = conda_activate
        self.max_cpus = int(min((cpu_count() / 4 * 3), 64))
        dir_anal = Path(analysis).resolve()

        # SOLAR YAML
        self.dict_inyml = yaml.safe_load(open(self.yaml_in))
        self.task_id = self.dict_inyml['task_id']
        self.dir_work = dir_anal.joinpath(str(self.dict_inyml['task_id']))

    def preparation(self):
        logging.info('准备步骤, 创建目录.')
        cml = f"mkdir -p {self.dir_work}/.rawdata {self.dir_work}/logs {self.dir_work}/Upload"
        logging.info("CommandLine: " + cml)
        run(cml, shell=True)

    def nextflow(self):
        logging.info('运行 Nextflow 流程.')
        cml = f"""
source {self.conda_activate} nextflow
cd {self.dir_work}
nextflow -log logs/nextflow.log run {self.dir_proj}/nf-wholegenome \
    -w ./.intermediate --threads {self.max_cpus} --csv .samplesheet.csv
"""
        logging.info("CommandLine: " + cml)
        run(cml, shell=True, executable='/bin/bash',
            capture_output=True, check=True)

    def upload(self):
        logging.info('复制结果到 Upload 目录.')
        for smp in self.dict_inyml['samples']:
            cml = f"""
cd {self.dir_work}/WholeGenomeAnalysis/{smp} && mkdir -p Upload/1.qc Upload/2.asm Upload/3.pred Upload/4.anno
# QC
cp -rf 1.qc/nanostat/{smp}.nanostat \
    1.qc/nanoplot/{smp}.*.png \
    1.qc/fastp \
    1.qc/fastqc \
    -t Upload/1.qc
# ASSEMBLY
cp -f 2.asm/unicycler/assembly.fasta Upload/{smp}.fa
cp -rf 2.asm/stat/{smp}_fa.stat.txt \
    2.asm/stat/{smp}.length.png \
    2.asm/depth/consistency.txt \
    2.asm/depth/consistency.all.txt \
    2.asm/depth/depth_base.stat.depth_GC.png \
    2.asm/checkm/checkm_parse.txt \
    2.asm/quast \
    -t Upload/2.asm
# PREDICTION
cp -f 3.pred/predict/predict.length.png \
    3.pred/predict/predict_kind.txt \
    3.pred/phispy/prophage_coordinates.txt \
    3.pred/island/island.txt \
    3.pred/repeat/INEs.tsv \
    3.pred/repeat/TandemRepeat.tsv \
    -t Upload/3.pred
# ANNOTATION
cp -f 4.anno/VFDB/Virulence.tsv \
    4.anno/VFDB/viru.fa \
    4.anno/VFDB/viru_ref.fa \
    4.anno/card/detail_card_resistance.txt \
    4.anno/card/detail_card_resistance.xlsx \
    4.anno/card/card.fa \
    4.anno/CAZy/CAZy_anno_stats.png \
    4.anno/CAZy/CAZY.txt \
    4.anno/eggmapper/COG/all.COG.bar.png \
    4.anno/eggmapper/GO/GO_anno_stats_level2.png \
    4.anno/eggmapper/GO/GO_anno_stats.xls \
    4.anno/eggmapper/KEGG/KEGG_anno_stats.png \
    4.anno/eggmapper/KEGG/KEGG_anno.stat.txt \
    4.anno/eggmapper/KEGG/KEGG_anno.txt \
    4.anno/swissprot/swissprot_result.tsv \
    4.anno/ResFinder/Resfinder_anno.txt \
    4.anno/ResFinder/ResFinder_Hit_in_genome_seq.fsa \
    4.anno/plasmid/plasmid.tsv \
    -t Upload/4.anno
# REPORT SOURCE
cp -rf {self.dir_proj}/assets/src Upload
# 240228 SOLAR 上传规范
mv {self.dir_work}/WholeGenomeAnalysis/{smp}/Upload {self.dir_work}/Upload/{smp}
"""
            logging.info("CommandLine: " + re.sub(" +", " ", cml))
            run(cml, shell=True)

    def report(self):
        logging.info('生成SOLAR报告.')
        for smp in self.dict_inyml['samples']:
            cml = f"perl {self.dir_proj}/bin/report_bacteria_denovo.pl {self.dir_work}/Upload/{smp} {smp}"
            logging.info("CommandLine: " + cml)
            run(cml, shell=True)

    def myzip(self):
        logging.info('压缩结果目录.')
        for smp in self.dict_inyml['samples']:
            cml = f"""
set -eu
cd {self.dir_work}/Upload
zip -r {smp}.zip {smp}
            """
            logging.info("CommandLine: " + cml)
            run(cml, shell=True)


class PipeV(Pipe):
    """ONT+ILMN混装病毒全基因组分析流程"""

    def __init__(self, inyaml, analysis, conda_activate) -> None:
        super().__init__(inyaml, analysis, conda_activate)

    def nextflow(self):
        logging.info('运行 Nextflow 流程.')
        cml = f"""
source {self.conda_activate} nextflow
cd {self.dir_work}
nextflow -log logs/nextflow.log run {self.dir_proj}/nf-wholegenome \
    -w ./.intermediate --threads {self.max_cpus} --biwf virassy --csv .samplesheet.csv
"""
        logging.info("CommandLine: " + cml)
        run(cml, shell=True, executable='/bin/bash',
            capture_output=True, check=True)

    def upload(self):
        logging.info('复制结果到 Upload 目录.')
        for smp in self.dict_inyml['samples']:
            cml = f"""
cd {self.dir_work}/WholeGenomeAnalysis/{smp} && mkdir -p Upload/1.qc Upload/2.asm Upload/3.pred Upload/4.anno
# QC
cp -rf 1.qc/nanostat/{smp}.nanostat \
    1.qc/nanoplot/{smp}.*.png \
    1.qc/fastp \
    1.qc/fastqc \
    -t Upload/1.qc
# ASSEMBLY
cp -f 2.asm/unicycler/assembly.fasta Upload/{smp}.fa
cp -rf 2.asm/stat/{smp}_fa.stat.txt \
    2.asm/stat/{smp}.length.png \
    2.asm/depth/consistency.txt \
    2.asm/depth/consistency.all.txt \
    2.asm/depth/depth_base.stat.depth_GC.png \
    2.asm/checkv/quality_summary.tsv \
    2.asm/quast \
    -t Upload/2.asm
# PREDICTION
cp -f 3.pred/predict/predict.length.png \
    3.pred/predict/predict_kind.txt \
    -t Upload/3.pred
# ANNOTATION
cp -f 4.anno/swissprot/swissprot_result.tsv \
    4.anno/KEGG/KEGG_anno.stat.tsv \
    4.anno/KEGG/KEGG_anno_stats.png \
    -t Upload/4.anno
# REPORT SOURCE
cp -rf {self.dir_proj}/assets/src Upload
# 240228 SOLAR Upload 规范
mv {self.dir_work}/WholeGenomeAnalysis/{smp}/Upload {self.dir_work}/Upload/{smp}
"""
            logging.info("CommandLine: " + re.sub(" +", " ", cml))
            run(cml, shell=True)

    def report(self):
        logging.info('生成SOLAR报告.')
        for smp in self.dict_inyml['samples']:
            cml = f"perl {self.dir_proj}/bin/report_virus_denovo.pl {self.dir_work}/Upload/{smp} {smp}"
            logging.info("CommandLine: " + cml)
            run(cml, shell=True)
