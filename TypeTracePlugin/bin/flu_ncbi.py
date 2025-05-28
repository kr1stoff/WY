#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/13 11:03
# @Last Modified by:   Ming
# @Last Modified time: 2022/6/13 11:03
from Bio import SeqIO
import logging
import click
import gzip

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

#### Some Global var
num2segment = {'1': "PB2", '2': "PB1", "3": "PA", "4": "HA",
               "5": "NP", "6": "NA", "7": "MP", "8": "NS"}


#### Some Function
def time_format(str_time):
    """
    Format the date
    """
    list_time = str_time.strip().split('/')
    time_length = len(list_time)
    if time_length == 1:
        res = f"{str_time}-01-01"
    elif time_length == 2:
        res = f"{list_time[0]}-{list_time[1]}-01"
    else:
        res = f"{list_time[0]}-{list_time[1]}-{list_time[2]}"
    return res


def get_isolate_name(list_line, _type):
    """

    """
    if '(' in list_line[7]:
        res = list_line[7].strip().split('(')[1].strip().split(')')[0]
    else:
        if list_line[1] == "Human":
            arr = [_type, list_line[4], '', list_line[5]]
        else:
            arr = [_type, list_line[1], list_line[4], '', list_line[5]]
        res = '/'.join(arr)
    res = res.replace(' ', '')
    return res


#### Main
@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
def cli():
    """
    NCBI 流感数据处理
    """
    pass


@click.command()
@click.option('-s', '--sequence',
              required=True,
              type=click.Path(),
              help="The influenza.fna.gz download from NCBI")
@click.option('--meta',
              required=True,
              type=click.Path(),
              help="The genome info file genomeset.dat.gz download from NCBI")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put fasta file")
def A(sequence, meta, out):
    """
    Prepare Influenza A download from NCBI https://ftp.ncbi.nih.gov/genomes/INFLUENZA/
    """
    logging.info(f"Parse the info file")
    info = {}
    with gzip.open(meta, 'rt') as IN:
        for line in IN:
            if 'Influenza A virus' in line:
                arr = line.strip().split('\t')
                isolate_name = get_isolate_name(arr, 'A')
                collection_date = time_format(arr[5])
                _type = f"A_/_{arr[3]}"
                # list_name = [f"{num2segment[arr[2]]}", isolate_name, collection_date, _type, arr[2]]
                list_name = [f"{num2segment[arr[2]]}", isolate_name]
                name = '|'.join(list_name)
                info[arr[0]] = name

    logging.info(f"Extract Sequence")
    with open(out, 'w') as OUT:
        for record in SeqIO.parse(gzip.open(sequence, 'rt'), 'fasta'):
            name = record.name.strip().split('|')[3]
            if name in info:
                print(f">{info[name]}|{name}\n{record.seq}", file=OUT)


@click.command()
@click.option('-s', '--sequence',
              required=True,
              type=click.Path(),
              help="The influenza.fna.gz download from NCBI")
@click.option('--meta',
              required=True,
              type=click.Path(),
              help="The genome info file genomeset.dat.gz download from NCBI")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put fasta file")
def B(sequence, meta, out):
    """
    Prepare Influenza B download from NCBI https://ftp.ncbi.nih.gov/genomes/INFLUENZA/
    """
    logging.info(f"Parse the info file")
    info = {}
    with gzip.open(meta, 'rt') as IN:
        for line in IN:
            if 'Influenza B virus' in line:
                arr = line.strip().split('\t')
                isolate_name = get_isolate_name(arr, 'B')
                collection_date = time_format(arr[5])
                _type = f"B_/_{arr[3]}"
                # list_name = [f"{num2segment[arr[2]]}", isolate_name, collection_date, _type, arr[2]]
                list_name = [f"{num2segment[arr[2]]}", isolate_name]
                name = '|'.join(list_name)
                info[arr[0]] = name

    logging.info(f"Extract Sequence")
    with open(out, 'w') as OUT:
        for record in SeqIO.parse(gzip.open(sequence, 'rt'), 'fasta'):
            name = record.name.strip().split('|')[3]
            if name in info:
                print(f">{info[name]}|{name}\n{record.seq}", file=OUT)


@click.command()
@click.option('-s', '--sequence',
              required=True,
              type=click.Path(),
              help="The influenza.fna.gz download from NCBI")
@click.option('--meta',
              required=True,
              type=click.Path(),
              help="The genome info file genomeset.dat.gz download from NCBI")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put fasta file")
def C(sequence, meta, out):
    """
    Prepare Influenza C download from NCBI https://ftp.ncbi.nih.gov/genomes/INFLUENZA/
    """
    logging.info(f"Parse the info file")
    info = {}
    with gzip.open(meta, 'rt') as IN:
        for line in IN:
            if 'Influenza C virus' in line:
                arr = line.strip().split('\t')
                isolate_name = get_isolate_name(arr, 'C')
                collection_date = time_format(arr[5])
                _type = f"C_/_{arr[3]}"
                # list_name = [f"{num2segment[arr[2]]}", isolate_name, collection_date, _type, arr[2]]
                list_name = [f"{num2segment[arr[2]]}", isolate_name]
                name = '|'.join(list_name)
                info[arr[0]] = name

    logging.info(f"Extract Sequence")
    with open(out, 'w') as OUT:
        for record in SeqIO.parse(gzip.open(sequence, 'rt'), 'fasta'):
            name = record.name.strip().split('|')[3]
            if name in info:
                print(f">{info[name]}|{name}\n{record.seq}", file=OUT)


cli.add_command(A)
cli.add_command(B)
cli.add_command(C)

if __name__ == "__main__":
    cli()
