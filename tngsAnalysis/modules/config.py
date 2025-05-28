from pathlib import Path
from collections import namedtuple
import yaml


def config():
    """读配置信息写到namedtuple中"""
    config_yaml = Path(__file__).parents[1].joinpath('conf/config.yaml')
    dict_conf = yaml.safe_load(open(config_yaml))
    Params = namedtuple("Params", dict_conf.keys())
    params = Params(**dict_conf)
    return params
