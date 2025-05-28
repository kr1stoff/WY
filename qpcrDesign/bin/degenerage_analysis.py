# @From: /sdbb/bioinfor/renchaobo/Develop/TracePlatform/bin/align_to_reference.py
import logging
from pathlib import Path
from subprocess import run
import click
import numpy as np
from Bio import SeqIO
import json
import pandas as pd
import pdb


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

def align(f_ref, f_seq, d_out, mafft):
    """
    多序列比对获取以参考为基准的对齐的序列
    """
    cmd = f"""set -e
{mafft} --quiet --6merpair --keeplength --addfragments {f_seq} {f_ref} > {d_out}/align.fa
"""
    run(cmd, shell=True)


def align_stat(info_snp, num_seq, f_align, f_stat, snp_ratio):
    """
    通过对齐的文件获取碱基突变信息

    :TODO 输出方式优化；添加阈值；直接输出兼并碱基
    """
    flag = 0
    with open(f_stat, 'w') as OUT:
        # {'A':[0,0,0,8,0,0,9,0,...], 
        #  'T':[0,9,0,0,0,7,0,0,...], 
        #   ......}
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
                            info_snp[alphabeta][i] -= 1 # [230915] 参考碱基的数量也要记录
                        else:
                            info_snp[alphabeta][i] += 1
            flag += 1
        # -------A------aA---------
        # ---t------T---T----------
        # .....
        # pdb.set_trace()
        for bs in info_snp.keys():
            obs = ''
            for idx in range(len(info_snp[bs])):
                # 0 : '-'
                if info_snp[bs][idx] == 0:
                    obs += '-'
                # 参考碱基那一条要算回去 -1 
                elif info_snp[bs][idx] < 0:
                    obs += '-'
                    info_snp[bs][idx] = abs(info_snp[bs][idx] - 1)
                # 小于snp阈值 : 小写字母
                elif info_snp[bs][idx] / num_seq < snp_ratio:
                    obs += bs.lower()
                # 大于snp阈值 : 大写字母
                else:
                    obs += bs
            print(obs, file=OUT)
    return info_snp


def write_base_number_json(info_snp:dict, outfile:str) -> None:
    """ 输出碱基频数矩阵"""
    fmt_info_snp = {k: list(v) for k,v in info_snp.items()}
    df_fmt_info_snp = pd.DataFrame(fmt_info_snp)
    df_fmt_info_snp.to_csv(outfile)


def get_degenerate_base(l_alpha):
    """
    转简并碱基, 突变不到阈值就小写参考碱基
    如果A碱基为参考碱基, 列表可能出现的情况
    A           :: 直接转简并
    A,T         :: 直接转简并
    A,T,c       :: 去掉小写转简并
    A,T,c,G     :: 去掉小写转简并
    A,t         :: 小写参考碱基
    A,t,c       :: 小写参考碱基
    """
    dgnrtdict = {('A',): 'A',
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
                    ('A', 'C', 'G', 'T'): 'N',
                    (): 'N', #参考碱基是N
                    }
    # 直接转简并
    if (len(l_alpha) == 1) or (''.join(l_alpha).isupper()):
        dgnrtbs = dgnrtdict[tuple(sorted(set(l_alpha)))]
    # 小写参考
    elif ''.join(l_alpha[1:]).islower():
        dgnrtbs = l_alpha[0].lower()
    # 去掉小写碱基后，转简并
    elif (not ''.join(l_alpha[1:]).islower()) and (not ''.join(l_alpha[1:]).isupper()):
        rmv_lwr = [bs for bs in l_alpha if bs.isupper()] 
        dgnrtbs = dgnrtdict[tuple(sorted(set(rmv_lwr)))]
    # 未知的组合
    else:
        logging.warning(f'未知组合: {l_alpha}')
    return dgnrtbs


def degenerate(f_stat, f_degenerate):
    """
    简并/小写碱基转化
    """
    with open(f_stat, 'r') as IN:
        ref = next(IN).strip()
        a = next(IN).strip()
        t = next(IN).strip()
        g = next(IN).strip()
        c = next(IN).strip()

    with open(f_degenerate, 'w') as OUT:
        print(ref, file=OUT)
        for i in range(len(ref)):
            # [230905] ref含有N碱基
            tmp = [bs for bs in [ref[i], a[i], t[i], g[i], c[i]] if bs not in ['-', 'N']]
            alpha = get_degenerate_base(list(tmp))
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
              default="/home/earthtest/miniconda3/envs/base_env/bin/mafft",
              help="The mafft path")
@click.option('--snp_ratio',
              default=0.05,
              type=float,
              show_default=True,
              help="The snp ratio to report")
def cli(sequence, name, out, mafft, snp_ratio):
    """
    简并碱基检查. 多序列比对, 按位置确定简并和SNP.
    """
    f_seq = Path(sequence).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True)

    d_prepare = d_out.joinpath("Prepare")
    d_prepare.mkdir(exist_ok=True)
    num_seq, len_ref = prepare(f_seq, name, d_prepare)
    info_snp = {"A": np.zeros(len_ref, dtype=int),
                "T": np.zeros(len_ref, dtype=int),
                "G": np.zeros(len_ref, dtype=int),
                "C": np.zeros(len_ref, dtype=int)}
    # 多序列比对
    d_align = d_out.joinpath("Align")
    d_align.mkdir(exist_ok=True)
    align(f"{d_prepare}/ref.fa", f"{d_prepare}/sequence.fa", d_align, mafft=mafft)
    # 统计简并
    f_stat = d_out.joinpath("result.txt")
    logging.debug(info_snp)
    align_stat(info_snp, num_seq, d_align.joinpath("align.fa"), f_stat, snp_ratio)
    # 碱基深度矩阵
    write_base_number_json(info_snp, d_out.joinpath('result.depth.csv'))
    # 简并标记
    f_degenerate = d_out.joinpath("degenerate.txt")
    degenerate(f_stat, f_degenerate)


if __name__ == "__main__":
    cli()
