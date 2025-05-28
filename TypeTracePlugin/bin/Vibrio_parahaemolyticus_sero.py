#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/11/1 11:16
# @Last Modified by:   Ming
# @Last Modified time: 2023/11/1 11:16
import logging
import re
from pathlib import Path
from subprocess import run

import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def run_blast(query: Path, database, out: Path, program: str,
              max_target_seqs: int = 1000,
              identity: float = 90,
              threads: int = 8,
              word_size: int = 15,
              **kwargs):
    """
    执行Blast比对分析

    :param query: The input query sequence
    :param database: The Vibrio parahaemolyticus database
    :param out: The output file
    :param program: The blastn path
    :param max_target_seqs: The number of result to show
    :param identity: The identity to use for blast
    :param threads: The thread number be used in blast
    :param word_size: Word size for wordfinder algorithm
    :param kwargs: blastn 所支持的参数
    """
    f_out = out.absolute()
    param = {"num_threads": threads,
             "word_size": word_size,
             "evalue": 1000,
             "dust": "no",
             "outfmt": "6",
             "perc_identity": identity,
             "max_target_seqs": max_target_seqs,
             "db": database,
             "query": query.absolute(),
             "out": f_out}

    for i, j in kwargs.items():
        param[i] = j

    logger.info(f"Blast against the database")
    cmd = ' '.join([program] + [f"-{k} {v}" for k, v in param.items()])
    logger.debug(cmd)
    run(cmd, shell=True, check=True)


def filter_blast(f_blast: Path, f_filter: Path, info_length, identity: float = 90, coverage: float = 90):
    """
    过滤blast的结果

    :param f_blast: The input blast m8 result
    :param f_filter: The output filter result
    :param info_length: The length info of database
    :param identity: The min identity
    :param coverage: The min coverage
    """
    with open(f_blast, 'r') as IN, open(f_filter, 'w') as OUT:
        for line in IN:
            arr = line.strip().split('\t')
            info = {"Query": arr[0], "Subject": arr[1], "Identity": float(arr[2]), "Alignment": int(arr[3]),
                    "Mismatch": int(arr[4]), "Gap": int(arr[5]), "QStart": int(arr[6]), "QEnd": int(arr[7]),
                    "SStart": int(arr[8]), "SEnd": int(arr[9]), "Evalue": float(arr[10]),
                    "BitScore": float(arr[11])}
            if (info["Identity"] >= identity) and (info["Alignment"] / info_length[info["Subject"]] * 100 >= coverage):
                print(*arr, sep="\t", file=OUT)


def _get_serotype(l_info: list):
    """
    获取单个的血清型
    """
    for i in l_info:
        if '_' in l_info:
            co_gene = re.sub('[ab]', lambda x: 'b' if x.group() == 'a' else 'a', i)
            if co_gene in set(l_info):
                name = i.strip().split('_')[0]
                return name
            else:
                continue
        else:
            return i


def predict_serotype(f_fiter: Path):
    """
    通过过滤后的blast结果预测血清型

    :param f_filter: The filtered blast result
    """
    info = {}
    with open(f_fiter, 'r') as IN:
        for line in IN:
            arr = line.strip().split("\t")
            key = tuple([arr[0], arr[1], arr[-1]])
            info[key] = int(arr[-1])

    pool = {"O_type": [], "K_type": []}
    for i, j in sorted(info.items(), key=lambda x: x[1], reverse=True):
        query, subject, score = i
        _type = subject.strip().split('*')[0]
        key = f"{_type[0]}_type"
        pool[key].append(_type)

    res = {"O_type": _get_serotype(pool["O_type"]),
           "K_type": _get_serotype(pool["K_type"])}
    return res


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-f', '--fasta',
              required=True,
              type=click.Path(),
              help="The input genome file of Vibrio parahaemolyticus")
@click.option('-d', '--database',
              required=True,
              type=click.Path(),
              help="The database dir for Vibrio parahaemolyticus")
@click.option('-l', '--length',
              required=True,
              type=click.Path(),
              help="The database seq length info file")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
@click.option('--blastn',
              default="/home/earthtest/miniconda3/envs/base_env/bin/blastn",
              show_default=True,
              type=click.Path(),
              help="The blastn path")
@click.option('--identity',
              default=90,
              show_default=True,
              type=float,
              help="The min identity")
@click.option('--coverage',
              default=90,
              show_default=True,
              type=float,
              help="The min coverage")
def main(fasta, database, length, out, blastn, identity, coverage):
    """
    副溶血弧菌血清型预测
    """
    f_fasta = Path(fasta).absolute()
    f_length = Path(length).absolute()
    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True, parents=True)

    logger.info(f"Run the blast against {database}")
    f_blast = d_out.joinpath("blast.m8")
    run_blast(f_fasta, database, f_blast, blastn)

    logger.info(f"Filter the blast result")
    info_length = {}
    with open(f_length, 'r') as IN:
        for line in IN:
            arr = line.strip().split("\t")
            info_length[arr[0]] = int(arr[1])
    f_filter = d_out.joinpath("filter.tsv")
    filter_blast(f_blast, f_filter, info_length, identity=identity, coverage=coverage)

    f_final = d_out.joinpath('result.tsv')
    logger.info(f"Get the serotype and output the result to {f_final}")
    # TODO: 还需要额外的判定条件，后期升级
    res = predict_serotype(f_filter)
    logger.debug(res)
    with open(f_final, 'w') as OUT:
        print(*["O_type", "K_type"], sep="\t", file=OUT)
        print(*[res["O_type"], res["K_type"]], sep="\t", file=OUT)


if __name__ == "__main__":
    main()
