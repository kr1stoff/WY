#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2023/2/20 9:32
# @Last Modified by:   Ming
# @Last Modified time: 2023/2/20 9:32
import logging
import click

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

#### Some Function





#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.0.0")
@click.option('-c', '--config',
              required=True,
              type=click.Path(),
              help="The YAML config file")
def main(config):
    """
    临时脚本
    """


    pass


if __name__ == "__main__":
    main()
