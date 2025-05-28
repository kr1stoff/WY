#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/14 20:17
# @Last Modified by:   Ming
# @Last Modified time: 2022/6/14 20:17
import logging
import os
from glob import glob
from pathlib import Path

import click
import yaml

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class Mutex(click.Option):
    """
    互斥参数
    """

    def __init__(self, *args, **kwargs):
        self.not_required_if: list = kwargs.pop("not_required_if")
        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (kwargs.get("help",
                                     "") + " [Option is mutually exclusive with " + ", ".join(
            self.not_required_if) + ".").strip() + "]"
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt: bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '" + str(
                            self.name) + "' is mutually exclusive with " + str(
                            mutex_opt) + ".")
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


class MyConfig(object):
    """
    配置文件对象
    """

    def __init__(self):
        """

        :param d_config: The project config dir
        """
        self._libpath = Path(__file__).parent
        self._configpath = self._libpath.parent.joinpath("config")
        self.config = {}
        for f_config in glob(str(self._configpath.joinpath("*yml"))):
            self.read(f_config)

    def read(self, fname, name=None):
        """
        配置文件的读取

        :param fname:
        :return:
        """
        if not name:
            name = Path(fname).stem
        logging.info(f"Add config {name} to object")
        self.config[name] = yaml.safe_load(open(fname, 'r').read())


class Flu(object):
    """
    流感扩增子流程
    """

    def __init__(self, sample, out, config, d_upload):
        """

        """
        self.d_out = os.path.abspath(out)
        os.makedirs(self.d_out, exist_ok=True)
        self.name = sample
        self.config = config.config
        self.cmd = ["set -e"]

        # 一些经常用到的变量
        self.bin = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin")
        logging.debug(self.bin)
        self.python = self.config['software']['python']
        self.perl = self.config['software']['perl']
        self.irma = os.path.join(os.path.dirname(self.bin), "IRMA/IRMA")
        self.rscript = self.config['software']['rscript']
        self.threads = self.config['project']['threads']
        self.memory = self.config['project']["memory"]
        self.db_all = self.config['database']["db_all"]
        self.db_ref = self.config['database']["db_ref"]
        self.segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

        # step
        self.qc()
        self.assemble()
        self.blast()
        self.phylotree()
        self.hostdetect()
        self.snp()
        self.coverage()
        self.report(d_upload)

    def timestamp(self, name, statu):
        """
        返回时间戳命令

        :param name: The Pipe name
        :return:
        """
        if statu == "start":
            cmd = f"echo -e [`date +%y-%m-%d.%H:%M:%S`] {name} Start"
        elif statu == "end":
            cmd = f"echo -e [`date +%y-%m-%d.%H:%M:%S`] {name} End\n"
        return cmd

    def qc(self):
        """
        质控

        TODO: 1. 添加绘图脚本
              2. 是否需要添加三代支持
              3. 搞成装饰器
        """
        self.d_qc = os.path.join(self.d_out, "01.QC")
        os.makedirs(self.d_qc, exist_ok=True)
        os.makedirs(f"{self.d_qc}/before", exist_ok=True)
        os.makedirs(f"{self.d_qc}/after", exist_ok=True)

        if self.config['project']["data_type"] == "SE" or self.config["project"]["data_type"] == "Nanopore":
            fq1 = self.config["project"]["samples"][self.name][0]
            if self.config["project"]["data_type"] == "Nanopore":
                self.config["project"]["data_type"] = "SE"
                self.config["project"]["filter"] == None
                f_fast5 = self.config["project"]["samples"][self.name][0]
                fq1 = f"{self.d_qc}/{self.name}_1.fq.gz"
                cmd = self.fast5(f_fast5, fq1)
                self.cmd.append(self.timestamp("Fast5 to Fastq", "Start"))
                self.cmd.append(cmd)
                self.cmd.append("Fast5 to Fastq", "End")
            cmd = self._se_qc(fq1, self.d_qc)
        elif self.config['project']["data_type"] == "PE":
            fq1 = self.config["project"]["samples"][self.name][0]
            fq2 = self.config["project"]["samples"][self.name][1]
            cmd = self._pe_qc(fq1, fq2, self.d_qc)
        self.cmd.append(self.timestamp("QC", "start"))
        self.cmd.append(cmd)
        cmd = f"{self.python} {self.bin}/fastp_json2table.py {self.d_qc}/{self.name}.json {self.d_qc}/qc.stat.tsv"
        self.cmd.append(cmd)
        self.cmd.append(self.timestamp("QC", "end"))

    def fast5(self, f_fast5, f_fastq):
        """
        Deal with Fast5 format to fastq from Nanopore
        """
        cmd = f"""unzip {f_fast5} -d {f_fast5}/tmp
{self.config['software']["guppy"]} --compress_fastq -c /opt/ont-guppy-cpu/data/dna_r9.4.1_450bps_hac.cfg -i {f_fast5}/tmp -s {f_fast5}/tmp2 --cpu_threads_per_caller 32 --num_callers 1
cat {f_fast5}/*fastq | gzip -c > {f_fastq}
"""
        return cmd

    def _se_qc(self, f_in, d_out):
        """
        SE QC

        :param f_in: The input fq1 read file
        :param d_out: The out put dir for the result
        :return: The cmd to exec
        """
        fastqc = self.config['software']['fastqc']
        fastp = self.config['software']["fastp"]
        bowtie2 = self.config['software']['bowtie2']
        cmd = f"""{fastqc} --quiet --extract -t {self.threads} -o {d_out}/before {f_in}
{fastp} -q 15 -u 40 -l 25 --thread 2 --cut_right --cut_window_size 20 --cut_mean_quality 30 -i {f_in} -o {d_out}/{self.name}_1.fq.gz -j {d_out}/{self.name}.json -h {d_out}/{self.name}.html
{fastqc} --quiet --extract -t {self.threads} -o {d_out}/after {d_out}/{self.name}_1.fq.gz"""
        if self.config["project"]["filter"]:
            filter_cmd = f"""\n{bowtie2} --local --mm -p {self.threads} -x {self.db_ref} -U {f_in} -S /dev/null --al {d_out}/{self.name}_1.filter.fq.gz 2>{d_out}/{self.name}.bowtie2.log"""
            cmd += filter_cmd
        return cmd

    def _pe_qc(self, fq1, fq2, d_out):
        """
        PE QC

        :param fq1: The input fq1 read file
        :param fq2: The input fq2 read file
        :param d_out: The out put dir for the result
        :return: The cmd to exec
        """
        fastqc = self.config['software']['fastqc']
        fastp = self.config['software']["fastp"]
        bowtie2 = self.config['software']['bowtie2']
        cmd = f"""{fastqc} --quiet --extract -t {self.threads} -o {d_out}/after {fq1} {fq2}
{fastp} -q 15 -u 40 -l 25 --thread 2 --cut_right --cut_window_size 20 --cut_mean_quality 30 -i {fq1} -I {fq2} -o {d_out}/{self.name}_1.fq.gz -O {d_out}/{self.name}_2.fq.gz -j {d_out}/{self.name}.json -h {d_out}/{self.name}.html
{fastqc} --quiet --extract -t {self.threads} -o {d_out}/after {d_out}/{self.name}_1.fq.gz {d_out}/{self.name}_2.fq.gz"""
        if self.config["project"]["filter"]:
            filter_cmd = f"""\n{bowtie2} --local --mm -p {self.threads} -x {self.db_ref} -1 {d_out}/{self.name}_1.fq.gz -2 {d_out}/{self.name}_2.fq.gz -S /dev/null --al-conc {d_out}/{self.name}_%.filter.fq.gz 2>{d_out}/{self.name}.bowtie2.log"""
            cmd += filter_cmd
        return cmd

    def assemble(self):
        """
        组装

        TODO: 1.暂时使用IRMA，待有时间搞清楚其原理后重写IRMA
        """
        self.d_assemble = os.path.join(self.d_out, "02.Assemble")
        os.makedirs(self.d_assemble, exist_ok=True)

        if self.config['project']["data_type"] == "SE":
            cmd = self._se_assemble(f"{self.d_assemble}/IRMA")
        elif self.config['project']["data_type"] == "PE":
            cmd = self._pe_assemble(f"{self.d_assemble}/IRMA")
        self.cmd.append(self.timestamp("Assemble", "start"))
        self.cmd.append(
            f"if [ -d {self.d_assemble}/IRMA ]; then rm -rf {self.d_assemble}/IRMA;fi")
        self.cmd.append(cmd)
        cmd = f"for i in `ls {self.d_assemble}/IRMA/*.fasta`;do name=$(basename $i .fasta | awk -F'_' '{{print $2}}'); ln -sf $i {self.d_assemble}/${{name}}.fasta;done"
        self.cmd.append(cmd)
        cmd = f"{self.python} {self.bin}/assemble_stat.py -d {self.d_assemble} -o {self.d_assemble}/assemble.stat.tsv"
        self.cmd.append(cmd)
        self.cmd.append(self.timestamp("Assemble", "end"))

    def _se_assemble(self, d_out):
        """
        SE 组装
        """
        fq = f"{self.d_qc}/{self.name}_1.filter.fq.gz" if self.config["project"][
            "filter"] else f"{self.d_qc}/{self.name}_1.fq.gz"
        cmd = f"bash {self.irma} FLU {fq} {d_out}"
        return cmd

    def _pe_assemble(self, d_out):
        """
        PE 组装
        """
        if self.config["project"]["filter"]:
            fq1 = f"{self.d_qc}/{self.name}_1.filter.fq.gz"
            fq2 = f"{self.d_qc}/{self.name}_2.filter.fq.gz"
        else:
            fq1 = f"{self.d_qc}/{self.name}_1.fq.gz"
            fq2 = f"{self.d_qc}/{self.name}_2.fq.gz"
        cmd = f"bash {self.irma} FLU {fq1} {fq2} {d_out}"
        return cmd

    def blast(self):
        """
        比对流感全库
        """
        self.d_blast = os.path.join(self.d_out, "03.Blast")
        os.makedirs(self.d_blast, exist_ok=True)
        self.cmd.append(self.timestamp("Blast", "start"))
        # param = f"-num_threads {self.threads} -evalue 1E-5 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp qcovus stitle' -max_target_seqs 50"
        param = f"-num_threads {self.threads} -evalue 1E-5 -outfmt 6 -max_target_seqs 50"
        cmd = f'{self.python} {self.bin}/flu_assemble_blast.py -i {self.d_assemble} -o {self.d_blast} --db {self.db_all} --blastn {self.config["software"]["blastn"]} --param "{param}"'
        self.cmd.append(cmd)
        self.cmd.append(self.timestamp("Blast", "end"))

    def phylotree(self):
        """
        提取序列并绘制进化树

        TODO: 更新整个流程所使用的Rscript路径
        """
        self.d_phylotree = os.path.join(self.d_out, "04.PhyloTree")
        os.makedirs(self.d_phylotree, exist_ok=True)
        self.cmd.append(self.timestamp("PhyloTree", "start"))
        cmd = f"""for f_blast in `ls {self.d_blast}/*blast.xls`
do
  ori_name=$(basename $f_blast .blast.xls)
  prefix={self.d_phylotree}/${{ori_name}}
  if [ -s $f_blast ]; then
    fullname=$(cut -f 1 $f_blast | head -n 1)
    cut -f 2 $f_blast | sort -u > ${{prefix}}.glist
    {self.perl} {self.bin}/fishInWinter.pl -ff fasta ${{prefix}}.glist {self.db_all}.fasta | cat - {self.d_assemble}/${{ori_name}}.fasta > ${{prefix}}.fasta
    {self.config['software']['mafft']} --thread {self.threads}  --maxiterate 1000 --quiet ${{prefix}}.fasta > ${{prefix}}.aln.fasta
    {self.config['software']['fasttree']} -quiet -nt ${{prefix}}.aln.fasta > ${{prefix}}.tre
    {self.rscript} {self.bin}/tree.R --branchlength --tree ${{prefix}}.tre --name ${{fullname}} --title ${{ori_name}} --prefix ${{prefix}}
  fi
done"""
        self.cmd.append(cmd)
        self.cmd.append(self.timestamp("PhyloTree", "end"))

    def hostdetect(self):
        """
        探究同近源序列宿主是否相同
        """
        self.d_hostdetect = os.path.join(self.d_out, "05.Hostdetect")
        os.makedirs(self.d_hostdetect, exist_ok=True)

    def snp(self):
        """
        突变分析
        """
        self.d_snp = os.path.join(self.d_out, "06.SNP")
        os.makedirs(self.d_snp, exist_ok=True)
        self.cmd.append(self.timestamp("SNP", "start"))
        cmd = f"""for i in {self.d_assemble}/IRMA/tables/*variants.txt
do
  ori_name=$(basename $i -variants.txt)
  name=$(basename $i -variants.txt | awk -F'_' '{{print $2}}')
  for j in variants insertions deletions
  do
    ln -sf  {self.d_assemble}/IRMA/tables/${{ori_name}}-${{j}}.txt {self.d_snp}/${{name}}.${{j}}.txt
  done
done"""
        self.cmd.append(cmd)
        self.cmd.append(self.timestamp("SNP", "end"))

    def coverage(self):
        """
        覆盖度分析

        TODO: 修改 genome_coverage.R 统计结果大于100%的bug
        """
        self.d_coverage = os.path.join(self.d_out, '07.Coverage')
        os.makedirs(self.d_coverage, exist_ok=True)
        samtools = self.config['software']['samtools']
        self.cmd.append(self.timestamp("Coverage", "start"))
        cmd = f"""for i in `ls {self.d_assemble}/IRMA/*.bam`
do
  if [ -s $i ]; then
    name=$(basename $i .bam | awk -F'_' '{{print $2}}')
    {samtools} sort -@ {self.threads} -o {self.d_coverage}/${{name}}.sort.bam $i
    {samtools} depth -a {self.d_coverage}/${{name}}.sort.bam > {self.d_coverage}/${{name}}.depth
    {self.rscript} {self.bin}/genome_coverage.R {self.d_coverage}/${{name}}.depth {self.d_coverage}/${{name}}
  fi
done
{self.python} {self.bin}/combain_coverage.py -d {self.d_coverage} -o {self.d_coverage}/coverage.tsv"""
        self.cmd.append(cmd)
        self.cmd.append(self.timestamp("Coverage", "end"))

    def report(self, d_report):
        """
        拷贝结果到指定目录并生成报告

        TODO: 05.Hostdetect

        :param d_report: The out put report dir
        :return cmd
        """
        self.d_report = os.path.abspath(d_report)
        os.makedirs(self.d_report, exist_ok=True)

        self.cmd.append(self.timestamp("Upload", "start"))
        cmd = f"""set +e
mkdir -p {self.d_report}/01.QC
if [ {self.config["project"]["data_type"]} = PE ]; then
  cp -rf {self.d_qc}/after/{self.name}_1*/Images/per_base_quality.png {self.d_report}/01.QC/qc_1.png
  cp -rf {self.d_qc}/after/{self.name}_2*/Images/per_base_quality.png {self.d_report}/01.QC/qc_2.png
else
  cp -rf {self.d_qc}/after/*/Images/per_base_quality.png {self.d_report}/01.QC/qc.png
fi
cp -rf {self.d_qc}/qc.stat.tsv {self.d_report}/01.QC
mkdir -p {self.d_report}/02.Assemble
cp -rf {self.d_assemble}/*fasta {self.d_report}/02.Assemble
cp -rf {self.d_assemble}/assemble.stat.tsv {self.d_report}/02.Assemble
mkdir -p {self.d_report}/03.Blast
cp -rf {self.d_blast}/*.annot.xls {self.d_report}/03.Blast
mkdir -p {self.d_report}/04.PhyloTree
cp -rf {self.d_phylotree}/*.tre {self.d_report}/04.PhyloTree
cp -rf {self.d_phylotree}/*.png {self.d_report}/04.PhyloTree
mkdir -p {self.d_report}/06.SNP
cp -rf {self.d_snp}/*.txt {self.d_report}/06.SNP
mkdir -p {self.d_report}/07.Coverage
cp -rf {self.d_coverage}/coverage.tsv {self.d_report}/07.Coverage
cp -rf {self.d_coverage}/*.genome_coverage_depth.png {self.d_report}/07.Coverage
cp -rf {self.d_coverage}/*.genome_coverage_depth_ylim1000.png {self.d_report}/07.Coverage
set -e"""
        self.cmd.append(cmd)
        self.cmd.append(self.timestamp("Upload", "end"))

        self.cmd.append(self.timestamp("Report and Package", "start"))
        self.cmd.append(f"{self.perl} {self.bin}/report.pl -outdir {self.d_report}")
        self.cmd.append(self.timestamp("Report and Package", "end"))

    def finish(self, f_script):
        """
        Finish the pipe and write command to script
        :return:
        """
        with open(f_script, 'w') as OUT:
            print('\n'.join(self.cmd), file=OUT)
