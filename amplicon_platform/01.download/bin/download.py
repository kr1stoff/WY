#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/12/6 14:52
# @Last Modified by:   Ming
# @Last Modified time: 2022/12/6 14:52
import logging
import subprocess
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
    genbank = genbank[~genbank["gbrs_paired_asm"].isin(refseq["# assembly_accession"])]
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
    index.update(dataframe[dataframe["taxid"] == taxonid].index)
    index.update(dataframe[dataframe["species_taxid"] == taxonid].index)
    if len(index) > 0:
        return dataframe.loc[list(index)]
    else:
        return None


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
@click.version_option(version="v1.0.0")
@click.option('--taxonid',
              required=True,
              type=click.Path(),
              help="The file contain taxonid to download")
@click.option('-o', "--out",
              required=True,
              type=click.Path(),
              help="The out put dir")
@click.option("--summary_refseq",
              required=True,
              type=click.Path(),
              help="The assembly_summary_refseq.txt")
@click.option("--summary_genbank",
              required=True,
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
def main(taxonid, out, summary_refseq, summary_genbank, thread, timeout, run, force, filter):
    """
    Download genome from NCBI's refseq and genbank database
    """
    filter_keyword = ["contaminated",
                      "chimeric",
                      "low quality sequence",
                      "misassembled",
                      "mixed culture",
                      "derived from environmental source",
                      "derived from metagenome",
                      "metagenome"]
    filter_keyword = '|'.join(filter_keyword)

    d_out = Path(out).absolute()
    d_out.mkdir(exist_ok=True)

    list_taxonid = set()
    logger.info(f"Parse the file: {taxonid}")
    with open(taxonid, 'r') as IN:
        for line in IN:
            line = int(line.strip())
            list_taxonid.add(line)

    logger.info(f"Read the summary file: {summary_refseq}")
    df_refseq = pd.read_csv(summary_refseq, sep="\t", skiprows=1, na_values=['na', ''])
    df_refseq = df_refseq[df_refseq["ftp_path"].notna()]
    if filter:
        df_refseq = df_refseq[~df_refseq["excluded_from_refseq"].str.contains(filter_keyword, na=False)]

    logger.info(f"Read the summary file: {summary_genbank}")
    df_genbank = pd.read_csv(summary_genbank, sep="\t", skiprows=1, na_values=['na', ''], low_memory=False)
    df_genbank = df_genbank[df_genbank["ftp_path"].notna()]
    if filter:
        df_genbank = df_genbank[~df_genbank["excluded_from_refseq"].str.contains(filter_keyword, na=False)]

    logger.info(f"Remove duplicate assem")
    df_summary = df_dup_remove(df_refseq, df_genbank)

    logger.info(f"Generate download script")
    f_script = d_out.joinpath("download.sh")
    l_download_script = []
    f_taxon_not_found = d_out.joinpath("not_found.list")
    l_taxon_not_found = []
    for i in list_taxonid:
        d_taxonid = d_out.joinpath(str(i))
        d_taxonid.mkdir(exist_ok=True)
        f_script_taxonid = d_taxonid.joinpath(f"{i}.sh")
        f_info_taxonid = d_taxonid.joinpath(f"{i}.info.tsv")
        df_query = query_summary(df_summary, i)
        if df_query is not None:
            df_query.to_csv(f_info_taxonid, sep="\t", index=False)
            try:
                cmd = df_query["ftp_path"].apply(get_download_cmd, d_out=d_taxonid)
            except:
                logger.error(df_query.loc[:, ["asm_name", "ftp_path"]])
                exit(-1)
            with open(f_script_taxonid, 'w') as OUT:
                for i in cmd:
                    if force:
                        print(i, file=OUT)
                    else:
                        if file_exist(i):
                            pass
                        else:
                            print(i, file=OUT)
            l_download_script.append(f_script_taxonid)
        else:
            l_taxon_not_found.append(i)
    with open(f_script, 'w') as OUT:
        for i in l_download_script:
            print(f"/usr/bin/parallel --timeout {timeout} -j {thread} < {i}", file=OUT)

    if len(l_taxon_not_found) > 0:
        logger.info(f"The taxonids not found are put into file: {f_taxon_not_found}")
        with open(f_taxon_not_found, 'w') as OUT:
            print(*l_taxon_not_found, sep="\n", file=OUT)

    if run:
        subprocess.run(f"bash {f_script}", shell=True)
    else:
        logger.info(f"Run the script {f_script} manual")


if __name__ == "__main__":
    main()
