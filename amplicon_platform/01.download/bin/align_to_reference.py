#!/use/bin/env python
# @Author: Chaobo Ren
# @Date:   2023-02-20 19:47:16
# @Last Modified by:   Chaobo Ren
# @Last Modified time: 2023-02-20 20:29:31
import logging
from pathlib import Path
from subprocess import run

import click
import numpy as np
from Bio import SeqIO

#### Some Functions ####
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def prepare(f_seq: Path, name: str, d_out: Path):
    """
    预处理
    """
    res = 0
    len_seq = 0
    f_ref = d_out.joinpath("ref.fa")
    f_other = d_out.joinpath("sequence.fa")
    pool = set()

    with open(f_ref, 'w') as OUT1, open(f_other, 'w') as OUT2:
        for record in SeqIO.parse(f_seq, "fasta"):
            if record.name not in pool:
                pool.add(record.name)
                res += 1
                if record.name == name:
                    len_seq = len(record.seq)
                    print(f">{record.name}\n{record.seq.upper()}", file=OUT1)
                else:
                    print(f">{record.name}\n{record.seq.upper()}", file=OUT2)
    return res, len_seq


def snp(f_ref, f_seq, d_out):
    """
    通过CALL变异来寻找兼并
    """
    cmd = f"""set -e
export PATH=/sdbb/bioinfor/renchaobo/Software/Samtools/samtools-1.15/bin:/sdbb/bioinfor/renchaobo/Software/minimap2/minimap2-2.24:/sdbb/bioinfor/renchaobo/Software/BCFtools/bcftools-1.15/bin:$PATH
/sdbb/bioinfor/renchaobo/Software/miniconda3/bin/python /sdbb/bioinfor/renchaobo/Software/MicroGMT/MicroGMT-v1.4/sequence_to_vcf.py -r {f_ref} -i assembly -fs {f_seq} -o {d_out} -kb
"""
    run(cmd, shell=True)


def align(f_ref, f_seq, d_out, mafft):
    """
    多序列比对获取以参考为基准的对齐的序列
    """
    cmd = f"""set -e
{mafft} --quiet --6merpair --keeplength --addfragments {f_seq} {f_ref} > {d_out}/align.fa
"""
    run(cmd, shell=True)


def align_stat(info_snp, number, f_align, f_stat, snp_ratio):
    """
    通过对齐的文件获取碱基突变信息

    :TODO 输出方式优化；添加阈值；直接输出兼并碱基
    """
    flag = 0
    with open(f_stat, 'w') as OUT:
        for record in SeqIO.parse(f_align, "fasta"):
            seq = record.seq.upper()
            if flag == 0:
                seq_ref = seq
                print(seq, file=OUT)
            else:
                for i in range(len(seq)):
                    alphabeta = seq[i]
                    if alphabeta in set(["A","T","G","C"]):
                        if seq[i] == seq_ref[i]:
                            pass
                        else:
                            info_snp[alphabeta][i] += 1
            flag += 1
        for i in info_snp.keys():
            # TODO：直接输出兼并碱基
            line = ['-' if j / number < snp_ratio else i for j in info_snp[i]]
            print(*line, sep="", file=OUT)


def stat(info_snp, number, f_ref, f_stat, d_snp: Path):
    """
    统计VCF信息
    """

    def deal_vcf(f_vcf):
        with open(f_vcf, 'r') as IN:
            for line in IN:
                if line.startswith('#'):
                    pass
                else:
                    arr = line.strip().split("\t")
                    pos = int(arr[1]) - 1
                    ori = arr[3]
                    mut = arr[4]
                    info_snp[mut][pos] += 1

    with open(f_stat, 'w') as OUT:
        for record in SeqIO.parse(f_ref, "fasta"):
            seq = record.seq
            print(f"{seq}", file=OUT)

        for f_vcf in d_snp.glob("*.vcf"):
            deal_vcf(f_vcf)

        for i in info_snp.keys():
            # TODO：直接输出兼并碱基
            # line = ['-' if j == 0 else f"{j / number:2f} * 100" for j in info_snp[i]]
            line = ['-' if j / number < 0.05 else i for j in info_snp[i]]
            print(*line, sep="", file=OUT)


dict_degenerate = {('A',): 'A',
                   ('C',): 'C',
                   ('G',): 'G',
                   ('T',): 'T',
                   ('A', 'G'): 'R',
                   ('C', 'T'): 'Y',
                   ('C', 'G'): 'S',
                   ('A', 'T'): 'W',
                   ('G', 'T'): 'K',
                   ('A', 'C'): 'M',
                   ('C', 'G', 'T'): 'B',
                   ('A', 'G', 'T'): 'D',
                   ('A', 'C', 'T'): 'H',
                   ('A', 'C', 'G'): 'V',
                   ('A', 'C', 'G', 'T'): 'N'}


def get_degenerate(l_alpha):
    return dict_degenerate[tuple(sorted(l_alpha))]


def degenerate(f_stat, f_degenerate):
    """

    """
    with open(f_stat, 'r') as IN:
        ref = next(IN).strip()
        a = next(IN).strip()
        t = next(IN).strip()
        g = next(IN).strip()
        c = next(IN).strip()
        # for line in IN:
        #     c = line.strip()

    with open(f_degenerate, 'w') as OUT:
        print(ref, file=OUT)
        for i in range(len(ref)):
            tmp = set([ref[i], a[i], t[i], g[i], c[i]])
            tmp.remove('-')
            alpha = get_degenerate(list(tmp))
            print(alpha, end="", file=OUT)


#####################
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-s', '--sequence',
              required=True,
              type=click.Path(),
              help="The input sequence file")
@click.option('-n', '--name',
              required=True,
              type=click.Path(),
              help="The reference name")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
@click.option('--mafft',
              required=True,
              type=click.Path(),
              show_default=True,
              default="/sdbb/bioinfor/zengql/tools/miniconda3/bin/mafft",
              help="The mafft path")
@click.option('--snp_ratio',
              default=0.05,
              type=float,
              show_default=True,
              help="The snp ratio to report")
def cli(sequence, name, out, mafft, snp_ratio):
    """
    获取突变信息，从而确定兼并碱基的设计
    """
    f_seq = Path(sequence).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True)

    d_prepare = d_out.joinpath("Prepare")
    d_prepare.mkdir(exist_ok=True)
    num_seq, len_ref = prepare(f_seq, name, d_prepare)
    # d_snp = d_out.joinpath("SNP")
    # d_snp.mkdir(exist_ok=True)
    info_snp = {"A": np.zeros(len_ref, dtype=int),
                "T": np.zeros(len_ref, dtype=int),
                "G": np.zeros(len_ref, dtype=int),
                "C": np.zeros(len_ref, dtype=int)}

    d_align = d_out.joinpath("Align")
    d_align.mkdir(exist_ok=True)
    align(f"{d_prepare}/ref.fa", f"{d_prepare}/sequence.fa", d_align, mafft=mafft)
    # snp(f"{d_prepare}/ref.fa", f"{d_prepare}/sequence.fa", d_snp)

    f_stat = d_out.joinpath("result.txt")
    logging.debug(info_snp)
    # stat(info_snp, num_seq, f"{d_prepare}/ref.fa", f_stat, d_snp)
    align_stat(info_snp, num_seq, d_align.joinpath("align.fa"), f_stat, snp_ratio)
    f_degenerate = d_out.joinpath("degenerate.txt")
    degenerate(f_stat, f_degenerate)


if __name__ == "__main__":
    cli()
