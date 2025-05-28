#!/usr/bin/env python


class Alignment():
    """新冠比对模块"""
    def __init__(self, _dict):
        self.dict_dirs = _dict['dict_dirs']
        self.prefix = _dict['prefix']
        self.fastq2 = _dict['fastq2']
        self.bwa = _dict['bwa']
        self.thread = _dict['thread']
        self.bwaidx = _dict['bwaidx']
        self.samtools = _dict['samtools']
        self.ivar = _dict['ivar']
        self.bed = _dict['bed']
        self.python = _dict['python']
        self.fgbio = _dict['fgbio']
        self.ktrim = _dict['ktrim']
        self.parallel_number = _dict['parallel_number']
        self.trim_software = _dict['trim_software']

    def bwa_samtools(self):
        """bwa + samtools 比对流程"""
        _fastq2 = f"{self.dict_dirs['qc']}/{self.prefix}.clean.2.fq" if self.fastq2 else " "
        cml = f"""
# 比对
{self.bwa} mem -t {self.thread} -M -Y -R '@RG\\tID:{self.prefix}\\tSM:{self.prefix}' \
    {self.bwaidx} {self.dict_dirs['qc']}/{self.prefix}.clean.1.fq {_fastq2} \
    | {self.samtools} view -@ 4 -hbS - \
    | {self.samtools} sort -@ 4 -o {self.dict_dirs['map']}/{self.prefix}.sorted.bam -
{self.samtools} index {self.dict_dirs['map']}/{self.prefix}.sorted.bam
        """
        return cml

    def trim_ivar(self):
        """ivar切引物. 三列bed, chrom start end"""
        cml = f"""
#~ 去引物
{self.ivar} trim -q 15 -m 30 -s 4 -e -b {self.bed} \
    -i {self.dict_dirs['map']}/{self.prefix}.sorted.bam \
    -p {self.dict_dirs['map']}/{self.prefix}.ivar_trim 
{self.samtools} sort -@ 4 \
    {self.dict_dirs['map']}/{self.prefix}.ivar_trim.bam \
    -o {self.dict_dirs['map']}/{self.prefix}.ivar_trim_sorted.bam
cp {self.dict_dirs['map']}/{self.prefix}.ivar_trim_sorted.bam {self.dict_dirs['map']}/{self.prefix}.bam
{self.samtools} index {self.dict_dirs['map']}/{self.prefix}.bam
            """
        return cml

    def trim_fgbio(self):
        """fgbio去引物. 五列table, chrom leftstart leftend rightstart rightend"""
        cml = f"""
# 去引物
{self.fgbio} TrimPrimers -p {self.bed} \
    -i {self.dict_dirs['map']}/{self.prefix}.sorted.bam \
    -o {self.dict_dirs['map']}/{self.prefix}.fgbio.bam
{self.samtools} sort -@ 4 \
    {self.dict_dirs['map']}/{self.prefix}.fgbio.bam \
    -o {self.dict_dirs['map']}/{self.prefix}.fgbio_sorted.bam
cp {self.dict_dirs['map']}/{self.prefix}.fgbio_sorted.bam {self.dict_dirs['map']}/{self.prefix}.bam
{self.samtools} index {self.dict_dirs['map']}/{self.prefix}.bam
        """
        return cml

    def myktrim(self):
        """自编去引物脚本, 五列table, chrom leftstart leftend rightstart rightend"""
        _fastq2 = f"-I {self.dict_dirs['qc']}/{self.prefix}.clean.2.fq" if self.fastq2 else " "
        cml = f"""
{self.python} {self.ktrim} -w {self.parallel_number} \
    -i {self.dict_dirs['qc']}/{self.prefix}.clean.1.fq {_fastq2} \
    -r {self.bwaidx} -p {self.bed} \
    -o {self.dict_dirs['qc']}/{self.prefix}.ktrim.fq
{self.bwa} mem -t {self.thread} -M -Y -R '@RG\\tID:{self.prefix}\\tSM:{self.prefix}' \
    {self.bwaidx} {self.dict_dirs['qc']}/{self.prefix}.ktrim.fq \
    | {self.samtools} view -@ 4 -hbS - \
    | {self.samtools} sort -@ 4 -o {self.dict_dirs['map']}/{self.prefix}.sorted.bam -
        """
        cml += self.copysorted_and_index()
        return cml

    def copysorted_and_index(self):
        """重复步骤, 复制.sorted.bam到.bam, 然后建索引"""
        cml = f"""
cp {self.dict_dirs['map']}/{self.prefix}.sorted.bam {self.dict_dirs['map']}/{self.prefix}.bam
{self.samtools} index {self.dict_dirs['map']}/{self.prefix}.bam
        """
        return cml
    
    def pipe(self):
        cml = ""
        if self.bed and self.trim_software == "ivar":
            cml += self.bwa_samtools()
            cml += self.trim_ivar()
        elif self.bed and self.trim_software == "fgbio":
            cml += self.bwa_samtools()
            cml += self.trim_fgbio()
        elif self.bed and self.trim_software == "ktrim":
            cml += self.myktrim()
        else:
            cml += self.bwa_samtools()
            cml += self.copysorted_and_index()
        return cml
