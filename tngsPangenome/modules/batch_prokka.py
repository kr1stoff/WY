import os
import sys
import logging
from pathlib import Path
from subprocess import run
sys.path.append(Path(__file__).parent)
from .config import config


params = config()


def batch_prokka(dir_fna, dir_gff, accessions_file):
    accessions = open(accessions_file).read().strip().split('\n')
    logging.debug(f'accessions length: {len(accessions)}')
    temp_batch_shell = Path(dir_gff).joinpath('temp_batch_prokka.sh')
    with open(temp_batch_shell, 'wt', newline='', encoding='utf-8') as g:
        for acc in accessions:
            if Path(f'{dir_gff}/{acc}.gff').exists(): continue #跑过了
            fna = list(Path(dir_fna).glob(f'{acc}*.fna'))[0]
            g.write(f'prokka --cpus {params.threads} --prefix {acc} --outdir {Path(dir_gff).joinpath(acc)} '\
                    f'--kingdom Bacteria --force --addgenes --quiet --locustag {acc} {fna} '\
                    f'&& mv {Path(dir_gff).joinpath(f"{acc}/{acc}.gff")} {Path(dir_gff).joinpath(f"{acc}.gff")} '\
                    f'&& rm -r {Path(dir_gff).joinpath(f"{acc}")}\n')
    cml = f'cat {temp_batch_shell} | parallel -j {params.prlnum}'
    logging.debug(cml)
    prokka_env = os.environ
    prokka_env['PATH'] = f'{Path(params.prokka).parent}:{prokka_env["PATH"]}'
    run(cml, shell=True, env=prokka_env)
    os.remove(temp_batch_shell)
    