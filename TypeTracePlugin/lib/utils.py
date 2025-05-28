#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/5/26 17:10
# @Last Modified by:   Ming
# @Last Modified time: 2022/5/26 17:10
import logging
import os
from pathlib import Path

import click
import yaml

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class Mutex(click.Option):
    """
    互斥参数
    """

    def __init__(self, *args, **kwargs):
        self.not_required_if: list = kwargs.pop("not_required_if")
        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (kwargs.get("help",
                                     "") + "Option is mutually exclusive with " + ", ".join(
            self.not_required_if) + ".").strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt: bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '" + str(
                            self.name) + "' is mutually exclusive with " + str(
                            mutex_opt) + ".")
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


class MyConfig(object):
    """
    一个包含了config文件夹内所有配置文件的对象
    """

    def __init__(self, ):
        """

        """
        self.config = {}

    def read(self, f_name, name=None):
        """
        从yaml文件中读取配置信息
        """
        f_name = Path(f_name)
        if not f_name.exists():
            logger.error(f"{f_name} not exists")
            raise FileNotFoundError
        if not name:
            name = f_name.stem
        self.config[name] = yaml.safe_load(open(f_name, 'r').read())

    def walk_through_config(self):
        """
        遍历config文件夹内的yaml文件
        """
        PATH = os.path.dirname(os.path.abspath(__file__))
