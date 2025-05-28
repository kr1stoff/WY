#!/use/bin/env python
# @Author: ChaoboRen
# @Date:   2022-01-18 18:52:28
# @Last Modified by:   ChaoboRen
# @Last Modified time: 2022-01-18 19:55:07
import logging
import os.path
import sys

import click

# My MetaBolome lib
Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, 'lib'))
from Pipe import Meta

__version__ = '1.0.0'
logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
LOG = logging.getLogger(__name__)

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-c', '--config',
              required=True,
              type=click.Path(),
              help="The yaml config file for pipe")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The out put dir")
@click.option('--run/--no-run',
              default="True",
              show_default=True,
              help="Whether direct run the batch")
def cli(config, out, run):
    """
    16S Pipe line for EARTH
    """
    config = os.path.abspath(config)
    out = os.path.join(out)

    logging.info(f"Start to build the 16S Pipe")
    project = Meta(config, out)
    project.prepare()
    project.qc()
    project.qiime2()
    project.taxon()
    project.diversity()
    project.beta()
    project.lefse()
    project.picrust2()
    project.package()
    project.finish(run=run)


if __name__ == "__main__":
    cli()
