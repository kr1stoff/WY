#!/usr/bin/env python

import re
import sys
import yaml
from pathlib import Path
from subprocess import run
import logging
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)


def muscle_seqtk_multi():
    """muscle多序列比对, seqkt排序, parallel并行"""
    fastas = Path(dir_fasta).glob('*.fa')
    hh = open(f'{dir_fasta}/all.sh', 'wt', encoding='utf-8', newline='')
    for fa in fastas:
        gh = open(f'{fa}.sh', 'wt', encoding='utf-8', newline='')
        gh.write(f"""
    {muscle} -in {fa} -out {fa}.1
    {seqkit} sort {fa}.1 | {seqkit} seq -w 0 > {fa}.3
        """)
        gh.close()
        hh.write(f'bash {fa}.sh\n')
    hh.close()

    cml = f"""
    cat {dir_fasta}/all.sh | {parallel} -j {parallel_num}
    """
    run(cml, shell=True)
    logging.debug(cml)

def make_core_gene():
    """orthologue 合并保守单拷贝基因序列"""
    seq_dict = {name:"" for name in dict_sample}
    for sof in Path(dir_fasta).glob('*.fa.3'):
        # 一样的单拷贝基因就不要了
        num_uniq_gene = len(set([str(record.seq) for record in SeqIO.parse(sof, "fasta")]))
        if num_uniq_gene == 1:
            continue
        for item in itertools.product(SeqIO.parse(sof, "fasta"), dict_sample.keys()):
            record, name = item
            if re.match(f"^{name}_", record.id):
                seq_dict[name] += str(record.seq)
    out_records = [SeqRecord(Seq(v), id=k, description=k) for k, v in seq_dict.items()]
    SeqIO.write(out_records, file_supergene, "fasta")

if __name__ == '__main__':
    # 输入
    if len(sys.argv) != 4:
        print('prog <dir_fasta> <inyaml> <file_supergene>')
        sys.exit(1)
    dir_fasta = sys.argv[1]
    inyaml = sys.argv[2]
    file_supergene = sys.argv[3]

    # 配置
    HOME = Path(__file__).parents[1]
    dict_soft = yaml.safe_load(open(f"{HOME}/conf/software.yml"))
    muscle = dict_soft['muscle']
    Gblocks = dict_soft['Gblocks']
    seqkit = dict_soft['seqkit']
    parallel = dict_soft['parallel']
    dict_params = yaml.safe_load(open(f"{HOME}/conf/parameters.yml"))
    parallel_num = dict_params['parallel']

    # 运行
    dict_input = yaml.safe_load(open(inyaml))
    dict_sample = dict_input['samples']
    muscle_seqtk_multi()
    make_core_gene()
