# -*- coding: utf-8 -*-
# @Author: Stone
# @Date:   2023-08-24 11:22

import logging
import os
import click

#### Some Functions ####
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
__version__ = '1.0.0'

########################
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    """
    Pipeline for result merge
    """
    pass


def get_result(file, col):
    result = {}
    with open(file, 'r') as IN:
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            per = arr[(int(col) - 1)]
            result[name] = per
    return result


# ===== *nested ===========================================================================================================

@click.command()
@click.option('-i', '--input',
              type=click.Path(),
              required=True,
              help='The input table')
@click.option('-o', '--output',
              type=click.Path(),
              default='./primer.merge.txt',
              help='The input table')
@click.option('-1', '--specificity_in',
              type=click.Path(),
              required=True,
              help='The specificity inprimer')
@click.option('-2', '--specificity_p',
              type=click.Path(),
              required=True,
              help='The specificity probe')
@click.option('-3', '--specificity_out',
              type=click.Path(),
              required=True,
              help='The specificity outprimer')
@click.option('-4', '--specificity_merge',
              type=click.Path(),
              required=True,
              help='The specificity merge')
@click.option('-5', '--inclusiveness_in',
              type=click.Path(),
              required=True,
              help='The inclusiveness inprimer')
@click.option('-6', '--inclusiveness_p',
              type=click.Path(),
              required=True,
              help='The inclusiveness probe')
@click.option('-7', '--inclusiveness_out',
              type=click.Path(),
              required=True,
              help='The inclusiveness outprimer')
@click.option('-8', '--inclusiveness_merge',
              type=click.Path(),
              required=True,
              help='The inclusiveness merge')
def nested(input, output, specificity_in, specificity_p, specificity_out, specificity_merge, inclusiveness_in, inclusiveness_p, inclusiveness_out, inclusiveness_merge):
    output = os.path.abspath(output)

    dic_specificity_in = get_result(specificity_in, 2)
    dic_specificity_p = get_result(specificity_p, 2)
    dic_specificity_out = get_result(specificity_out, 2)
    dic_specificity_merge = get_result(specificity_merge, 2)

    dic_inclusiveness_in = get_result(inclusiveness_in, 2)
    dic_inclusiveness_p = get_result(inclusiveness_p, 2)
    dic_inclusiveness_out = get_result(inclusiveness_out, 2)
    dic_inclusiveness_merge = get_result(inclusiveness_merge, 4)
    with open(input, 'r') as IN, open(output, 'w') as OUT:
        print(*['ID', 'In_specificity', 'Probe_specificity', 'Out_specificity', 'In_inclusiveness', 'Probe_inclusiveness',
              'Out_inclusiveness', 'Merge_specificity', 'Merge_inclusiveness'], sep='\t', file=OUT)
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            dic_specificity_in.setdefault(name, 0)
            dic_specificity_p.setdefault(name, 0)
            dic_specificity_out.setdefault(name, 0)
            dic_inclusiveness_in.setdefault(name, 0)
            dic_inclusiveness_p.setdefault(name, 0)
            dic_inclusiveness_out.setdefault(name, 0)
            dic_specificity_merge.setdefault(name, 0)
            dic_inclusiveness_merge.setdefault(name, 0)
            print(*[name, dic_specificity_in[name], dic_specificity_p[name], dic_specificity_out[name], dic_inclusiveness_in[name],
                  dic_inclusiveness_p[name], dic_inclusiveness_out[name], dic_specificity_merge[name], dic_inclusiveness_merge[name]], sep='\t', file=OUT)


# ===== *qPCR ===========================================================================================================

@click.command()
@click.option('-i', '--input',
              type=click.Path(),
              required=True,
              help='The input table')
@click.option('-o', '--output',
              type=click.Path(),
              default='./primer.merge.txt',
              help='The input table')
@click.option('-1', '--specificity_in',
              type=click.Path(),
              required=True,
              help='The specificity inprimer')
@click.option('-2', '--specificity_p',
              type=click.Path(),
              required=True,
              help='The specificity probe')
@click.option('-3', '--inclusiveness_in',
              type=click.Path(),
              required=True,
              help='The inclusiveness inprimer')
@click.option('-4', '--inclusiveness_p',
              type=click.Path(),
              required=True,
              help='The inclusiveness probe')
def qpcr(input, output, specificity_in, specificity_p, inclusiveness_in, inclusiveness_p):
    output = os.path.abspath(output)

    dic_specificity_in = get_result(specificity_in, 2)
    dic_specificity_p = get_result(specificity_p, 2)

    dic_inclusiveness_in = get_result(inclusiveness_in, 2)
    dic_inclusiveness_p = get_result(inclusiveness_p, 2)
    with open(input, 'r') as IN, open(output, 'w') as OUT:
        print(*['ID', 'In_specificity', 'Probe_specificity',
              'In_inclusiveness', 'Probe_inclusiveness'], sep='\t', file=OUT)
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            dic_specificity_in.setdefault(name, 0)
            dic_specificity_p.setdefault(name, 0)
            dic_inclusiveness_in.setdefault(name, 0)
            dic_inclusiveness_p.setdefault(name, 0)
            print(*[name, dic_specificity_in[name], dic_specificity_p[name],
                  dic_inclusiveness_in[name], dic_inclusiveness_p[name]], sep='\t', file=OUT)

# ===== *tNGS ===========================================================================================================


@click.command()
@click.option('-i', '--input',
              type=click.Path(),
              required=True,
              help='The input table')
@click.option('-o', '--output',
              type=click.Path(),
              default='./primer.merge.txt',
              help='The input table')
@click.option('-1', '--specificity_out',
              type=click.Path(),
              required=True,
              help='The specificity outprimer')
@click.option('-2', '--inclusiveness_out',
              type=click.Path(),
              required=True,
              help='The inclusiveness outprimer')
def tngs(input, output, specificity_out, inclusiveness_out):
    output = os.path.abspath(output)

    dic_specificity_out = get_result(specificity_out, 2)

    dic_inclusiveness_out = get_result(inclusiveness_out, 2)
    with open(input, 'r') as IN, open(output, 'w') as OUT:
        print(*['ID', 'Out_specificity', 'Out_inclusiveness'], sep='\t', file=OUT)
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            dic_specificity_out.setdefault(name, 0)
            dic_inclusiveness_out.setdefault(name, 0)
            print(*[name, dic_specificity_out[name],
                  dic_inclusiveness_out[name]], sep='\t', file=OUT)


cli.add_command(nested)
cli.add_command(qpcr)
cli.add_command(tngs)


if __name__ == "__main__":
    cli()
