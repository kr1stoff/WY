#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/6/14 20:01
# @Last Modified by:   Ming
# @Last Modified time: 2022/6/14 20:01
import logging
import subprocess
import sys
from pathlib import Path

import click
import yaml

#### My own lib and some global setting
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
LIBPATH = Path(__file__).parent.joinpath("lib")
BIN = Path(__file__).parent.joinpath("bin")
sys.path.append(str(LIBPATH))
from Pipe import Mutex
from Pipe import MyConfig
from Pipe import Flu


#### Some Function
def venus2yaml(yaml_venus, yaml_pipe):
    """
    将VENUS生成的配置文件转为流程所需配置文件

    :param yaml_venus: The yaml file generate by VENUS
    :param yaml_pipe: The yaml file needed by the pipe
    """
    cfg_venus = yaml.safe_load(open(yaml_venus, 'r').read())
    res = {}
    res["project"] = cfg_venus["task_id"]
    res["db_all"] = False
    res["db_ref"] = False
    res["primer"] = cfg_venus["primer_seq"]
    res["filter"] = False
    if cfg_venus["pe_or_se"] == "单端":
        res["data_type"] = "SE"
    elif cfg_venus["pe_or_se"] == "双端":
        res["data_type"] = "SE"
    elif cfg_venus["pe_or_se"] == "Nanopore":
        res["data_type"] = "Nanopore"
    else:
        raise TypeError(f"不支持的数据类型")
    res["samples"] = {}
    for sample in cfg_venus["samples"].keys():
        if res["data_type"] == "SE":
            res["samples"][sample] = [cfg_venus["samples"][sample]["sequencing_data1"]]
        else:
            res["samples"][sample] = [cfg_venus["samples"][sample]["sequencing_data1"],
                                      cfg_venus["samples"][sample]["sequencing_data2"]]

    res["threads"] = 30
    res["memory"] = 100
    res["parallel"] = 8
    with open(yaml_pipe, 'w') as OUT:
        yaml.dump(res, OUT, default_flow_style=False)


#### Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="v1.1.4")
@click.option('--venus',
              cls=Mutex,
              not_required_if=["config", "out"],
              type=click.Path(),
              help="The yaml config file generate by VENUS")
@click.option('-c', '--config',
              cls=Mutex,
              not_required_if=["venus"],
              type=click.Path(),
              help="The YAML config file for Flu amplicon")
@click.option('-o', '--out',
              cls=Mutex,
              not_required_if=["venus"],
              type=click.Path(),
              help="The out put dir")
@click.option('--run/--no-run',
              default=True,
              show_default=True,
              type=click.BOOL,
              help="Whether run the pipeline directly")
def main(venus, config, out, run):
    """
    流感扩增子流程
    """
    logging.info(f"Start VENUS Influenza Amplicon Pipeline")
    if venus:
        f_venus = Path(venus).absolute()
        f_config = f_venus.parent.joinpath("config.yml")
        venus2yaml(f_venus, f_config)
        d_out = f_venus.parent
    else:
        f_config = Path(config)
        d_out = Path(out)
    d_out.mkdir(exist_ok=True)
    d_upload = d_out.joinpath("Upload")
    d_upload.mkdir(exist_ok=True)
    obj_config = MyConfig()
    obj_config.read(f_config, name="project")

    d_script = d_out.joinpath("Shell")
    d_script.mkdir(exist_ok=True)
    list_sample_script = []
    for sample in obj_config.config["project"]["samples"].keys():
        d_sample = str(d_out.joinpath(sample).absolute())
        d_sample_upload = d_upload.joinpath(sample)
        d_sample_upload.mkdir(exist_ok=True)
        d_sample_script = d_script.joinpath(sample)
        d_sample_script.mkdir(exist_ok=True)
        f_script = d_sample_script.joinpath(f"{sample}.sh")
        pipe = Flu(sample, d_sample, obj_config, d_sample_upload)
        pipe.finish(f_script)
        list_sample_script.append(f_script)

    f_all_script = d_script.joinpath("all.sh")
    with open(f_all_script, 'w') as OUT:
        for i in list_sample_script:
            print(f"bash {str(i.absolute())} > {str(i.absolute())}.log 2>&1", file=OUT)

    if run:
        cmd = f"{obj_config.config['software']['parallel']} -j {obj_config.config['project']['parallel']} < {f_all_script.absolute()} > {f_all_script.absolute()}.log 2>&1;echo 0"
        logging.debug(cmd)
        subprocess.run(cmd, shell=True)


if __name__ == "__main__":
    main()
