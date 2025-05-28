#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/12/6 14:52
# @Last Modified by:   Ming
# @Last Modified time: 2023/08/24 14:52
import logging
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse

import click
import pandas as pd

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


#### Some Function
def df_dup_remove(refseq, genbank):
    """

    :param refseq: The dataframe of refseq summary
    :param genbank: The dataframe of genbank summary
    :return:
    """
    genbank = genbank[~genbank["gbrs_paired_asm"].isin(refseq["assembly_accession"])]
    res = pd.concat([refseq, genbank]).reset_index(drop=True)
    return res


def query_summary(dataframe, taxonid):
    """
    Get the download path of the given taxon id

    :param dataframe: The Dataframe of summary info
    :param taxonid: The query taxon id

    :return: A list of the download path
    """
    res = []
    index = set()
    for i in taxonid:
        index.update(dataframe[dataframe["taxid"] == int(i)].index)
        index.update(dataframe[dataframe["species_taxid"] == int(i)].index)
    if len(index) > 0:
        return dataframe.loc[list(index)]
    else:
        return None


def ascp_download(url, d_out: Path,
                  ascp="/home/renchaobo/.aspera/connect/bin/ascp",
                  key="/home/renchaobo/.aspera/connect/etc/asperaweb_id_dsa.openssh"):
    """
    使用ascp下载文件
    """
    ori_path = urlparse(url)
    stem_path = ori_path.path
    name = stem_path.strip().split('/')[-1]
    cmd = f"{ascp} -i {key} -k 1 -T -l100m anonftp@ftp.ncbi.nlm.nih.gov:{url} {d_out}"
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode == 0:
        return d_out.joinpath(name)
    else:
        sys.exit(logger.error(f"{ret.returncode}: {cmd}"))


def expand_taxon(taxon,
                 taxonkit="/home/renchaobo/.local/bin/taxonkit",
                 d_data="/home/renchaobo/.taxonkit"):
    """
    将taxonid展开，获取其子id
    """
    res = set()
    cmd = f"{taxonkit} list --data-dir {d_data} -i {taxon}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    for i in ret.stdout.decode().strip().split("\n"):
        if i.strip():
            res.add(i.strip())
    return res


def get_download_cmd(ori_path, d_out,
                     ascp="/home/renchaobo/.aspera/connect/bin/ascp",
                     key="/home/renchaobo/.aspera/connect/etc/asperaweb_id_dsa.openssh"):
    """
    Get the download path from the ftp path

    :param oripath: The ori path in summary file
    :param d_out: The output dir
    :param ascp: The ascp path
    :param key: The ascp key

    :return: The ascp path to download
    """
    ori_path = urlparse(ori_path)
    stem_path = ori_path.path
    name = stem_path.strip().split('/')[-1]
    cmd = f"{ascp} -i {key} -k 1 -T -l100m anonftp@ftp.ncbi.nlm.nih.gov:{stem_path}/{name}_genomic.fna.gz {d_out}"
    return cmd


def file_exist(cmd):
    """
    过滤已下载完成的基因组的下载命令

    :param cmd: The download cmd
    """
    download_path, d_out = cmd.strip().split(':')[-1].strip().split()
    f_name = download_path.strip().split('/')[-1]
    f_out = Path(f"{d_out}/{f_name}")
    if f_out.exists():
        return True
    else:
        return False


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--taxonid',
              required=True,
              type=click.Path(),
              help="The taxonid's sep by ',' to download")
@click.option('--ptype',
              required=True,
              type=click.Choice(['Bacteria', "Viruses", "Fungi"]),
              help="The pathogen type")
@click.option('--db',
              required=True,
              type=click.Path(),
              help="The database dir")
@click.option('--name',
              required=True,
              type=click.Path(),
              help="The dir name for species database to store")
@click.option('--latin',
              type=click.Path(),
              required=True,
              help="The latin name")
@click.option("--summary_refseq",
              type=click.Path(),
              help="The assembly_summary_refseq.txt")
@click.option("--summary_genbank",
              type=click.Path(),
              help="The assembly_summary_genbank.txt")
@click.option("--thread",
              default=24,
              show_default=True,
              type=int,
              help="The parallel number to download the genomes")
@click.option("--timeout",
              default=45,
              type=int,
              show_default=True,
              help="The time out for parallel when download genome")
@click.option("--run/--no-run",
              default=False,
              show_default=True,
              help="Whether run the download script directly")
@click.option("--force/--no-force",
              default=False,
              show_default=True,
              help="Whether re download the genome even it is already exists")
@click.option("--filter/--no-filter",
              default=True,
              show_default=True,
              help="Whether filter the genome by column excluded_from_refseq in summary file")
@click.option("--ascp",
              default="/home/renchaobo/.aspera/connect/bin/ascp",
              show_default=True,
              help="The ascp path")
@click.option("--ascpkey",
              default="/home/renchaobo/.aspera/connect/etc/asperaweb_id_dsa.openssh",
              show_default=True,
              help="The ascp key")
@click.option("--taxonkit",
              default="/home/renchaobo/.local/bin/taxonkit",
              show_default=True,
              help="The taxonkit path")
@click.option("--taxondb",
              default="/home/renchaobo/.taxonkit",
              show_default=True,
              help="The taxonkit database path")
@click.option("--parallel",
              default="/usr/bin/parallel",
              show_default=True,
              help="The parallel path")
@click.option("--python",
              default="/sdbb/bioinfor/renchaobo/Software/miniconda3/bin/python",
              show_default=True,
              help="The python path")
def main(taxonid, ptype, db, name, latin, summary_refseq, summary_genbank, thread, timeout, run, force, filter,
         ascp, ascpkey, taxonkit, taxondb, parallel, python):
    """
    Download genome from NCBI's refseq and genbank database and Nucleotide
    """
    d_database = Path(db).absolute()
    d_database.mkdir(exist_ok=True, parents=True)

    # assembly_summary_refseq.txt 文件处理
    if summary_refseq:
        f_refseq = Path(summary_refseq).absolute()
    elif d_database.joinpath("assembly_summary_refseq.txt").exists():
        f_refseq = d_database.joinpath("assembly_summary_refseq.txt")
    else:
        logger.info(f"Download assembly_summary_refseq.txt")
        f_refseq = ascp_download("/genomes/refseq/assembly_summary_refseq.txt", d_database, ascp=ascp, key=ascpkey)

    # assembly_summary_genbank.txt 文件处理
    if summary_genbank:
        f_genbank = Path(summary_genbank).absolute()
    elif d_database.joinpath("assembly_summary_genbank.txt").exists():
        f_genbank = d_database.joinpath("assembly_summary_genbank.txt")
    else:
        logger.info(f"Download assembly_summary_genbank.txt")
        f_genbank = ascp_download("/genomes/genbank/assembly_summary_genbank.txt", d_database, ascp=ascp, key=ascpkey)

    # 数据下载
    d_out = d_database.joinpath(name)
    d_out.mkdir(exist_ok=True, parents=True)

    l_taxonid = set()
    logger.info(f"Expand TaxonID {taxonid}")
    for i in taxonid.strip().split(','):
        l_taxonid.update(expand_taxon(i, taxonkit=taxonkit, d_data=taxondb))

    filter_keyword = ["contaminated",
                      "chimeric",
                      "low quality sequence",
                      "misassembled",
                      "mixed culture",
                      "derived from environmental source",
                      "derived from metagenome",
                      "metagenome"]
    filter_keyword = '|'.join(filter_keyword)

    logger.info(f"Read the summary file: {f_refseq}")
    df_refseq = pd.read_csv(f_refseq, sep="\t", skiprows=1, na_values=['na', ''], low_memory=False)
    df_refseq.rename(columns={df_refseq.columns[0]: "assembly_accession"}, inplace=True)
    df_refseq = df_refseq[df_refseq["ftp_path"].notna()]
    if filter:
        df_refseq = df_refseq[~df_refseq["excluded_from_refseq"].str.contains(filter_keyword, na=False)]

    logger.info(f"Read the summary file: {f_genbank}")
    df_genbank = pd.read_csv(f_genbank, sep="\t", skiprows=1, na_values=['na', ''], low_memory=False)
    df_genbank.rename(columns={df_genbank.columns[0]: "assembly_accession"}, inplace=True)
    df_genbank = df_genbank[df_genbank["ftp_path"].notna()]
    if filter:
        df_genbank = df_genbank[~df_genbank["excluded_from_refseq"].str.contains(filter_keyword, na=False)]

    logger.info(f"Remove duplicate assem")
    df_summary = df_dup_remove(df_refseq, df_genbank)

    logger.info(f"Generate download script")
    d_seq = d_out.joinpath("seq")
    d_seq.mkdir(exist_ok=True, parents=True)
    f_download = d_seq.joinpath("download.sh")
    f_metadata = d_out.joinpath("metadata.tsv")
    df_query = query_summary(df_summary, l_taxonid)
    if df_query is not None:
        df_query.to_csv(f_metadata, sep="\t", index=False)
    try:
        cmd = df_query["ftp_path"].apply(get_download_cmd, d_out=d_seq)
    except:
        sys.exit(logger.error(df_query.loc[:, ["asm_name", "ftp_path"]]))
    with open(f_download, 'w') as OUT:
        for i in cmd:
            if force:
                print(i, file=OUT)
            else:
                if file_exist(i):
                    pass
                else:
                    print(i, file=OUT)

    f_script = d_out.joinpath("download.sh")
    with open(f_script, 'w') as OUT:
        print(f"{parallel} --timeout {timeout} -j {thread} < {f_download}", file=OUT)
        if ptype == "Viruses":
            d_bin = Path(__file__).absolute().parent
            print(
                f"{python} {d_bin}/entrez_search.py '({latin})' '(complete genome)' -o {d_seq}/list.txt", file=OUT)
            print(
                f"{python} {d_bin}/entrez_download.py -t {d_seq}/list.txt -p {d_seq}/nucleotide", file=OUT)

    if run:
        logger.info(f"Start to download...")
        subprocess.run(f"bash {f_script} > {f_script}.log 2>&1", shell=True)
        logger.info(f"Finish download")
    else:
        logger.info(f"Run the script {f_script} manual")


if __name__ == "__main__":
    main()
