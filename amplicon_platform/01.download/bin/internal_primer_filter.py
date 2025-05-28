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

__version__ = '1.0.0'

import yaml

logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-f', '--file',
              required=True,
              type=click.Path(),
              help="The Primer3 out put internal primer result")
@click.option('-o', '--out',
              required=True,
              type=click.Path(),
              help="The output dir")
def main(file, out):
    """
    内引物结果过滤
    """
    log
