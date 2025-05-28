#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/9 9:24
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/9 9:24
import logging
import os.path
import sys
from pathlib import Path

import click

__version__ = '1.1.0'

import yaml

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())
logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option("--run/--no-run",
              default=True,
              show_default=True,
              help="Whether run the pipeline or just generate the script")
@click.pass_context
def cli(ctx, run):
    """
    微远公卫引物设计平台
    """
    ctx.ensure_object(dict)
    ctx.obj["run"] = run


@click.command(context_settings=CONTEXT_SETTINGS)
@click.pass_context
@click.option('-c', '--config',
              required=True,
              type=click.Path(),
              help="The config for do multi PCR primer design")
@click.option('-o', '--out',
              required=True,
              help="The out put dir")
def mpcr(ctx, config, out):
    """
    多重PCR引物设计流程
    """
    logger.info(f"Load the config file {config}")
    f_config = os.path.abspath(config)
    config = yaml.safe_load(open(f_config))

    from lib.mPCR import Primer
    logger.info("Start the multi PCR primer design program")
    reference = os.path.abspath(config["reference"])
    genomes = os.path.abspath(config["genomes"])
    d_out = os.path.abspath(out)
    project = Primer(reference, genomes, d_out, run=ctx.obj["run"])
    project.format()
    project.snp()
    project.mask()
    project.design(f_config=f_config)
    project.evaluate()
    project.finish()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.pass_context
@click.option('-c', '--config',
              required=True,
              type=click.Path(),
              help="The config for pre-amplification primer design")
@click.option('-o', '--out',
              required=True,
              help="The out put dir")
def ppcr(ctx, config, out):
    """
    预扩增平台引物设计
    """
    f_config = Path(config).absolute()
    logger.info(f"Load the config file {f_config}")
    config = yaml.safe_load(open(f_config))

    from lib.pPCR import PPrimer
    logger.info("Start pre-amplification primer design")
    reference = Path(config["reference"]).absolute()
    genomes = Path(config["genomes"]).absolute()
    d_out = Path(out).absolute()
    project = PPrimer(reference, genomes, d_out, run=ctx.obj["run"])
    project.format()
    project.snp(freq=float(config["snp_freq"]))
    project.mask()
    project.internal_design(f_config=f_config)
    project.external_design()
    project.merge()
    project.finish()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.pass_context
@click.option('-c', '--config',
              required=True,
              type=click.Path(),
              help="The config for sequence primer design")
@click.option('-o', '--out',
              required=True,
              help="The out put dir")
def sequence(ctx, config, out):
    """
    测序引物设计流程
    """
    f_config = Path(config).absolute()
    d_out = Path(out).absolute()
    logger.info(f"Load the config file {f_config}")
    config = yaml.safe_load(open(f_config))

    from lib.Sequence import Sequence
    logger.info("Start Sequence Primer design")
    reference = Path(config["reference"]).absolute()
    genomes = Path(config["genomes"]).absolute()
    project = Sequence(reference, genomes, d_out, run=ctx.obj["run"])
    project.prepare()
    project.snp(freq=float(config["freq"]))
    project.design(Path(config["region"]).absolute(), f_config, freq=float(config["freq"]))
    project.specificity(config["taxonid"])
    project.inclusive(config["inclusive_database"])
    project.merge(config["taxonname"])
    project.finish()


cli.add_command(mpcr)
cli.add_command(ppcr)
cli.add_command(sequence)

if __name__ == "__main__":
    cli()
