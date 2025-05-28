#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/9 9:26
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/9 9:26
import logging
import click

__version__ = '1.0.0'
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
@click.option('-c', '--config',
              required=True,
              type=click.Path(),
              help="The YAML config file for metabolome")
def main(config):
    """
    DESC
    """
    pass


if __name__ == "__main__":
    main()
