#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Ming
# @Date:   2019-12-08 23:36:30
# @Last Modified by:   Ming
# @Last Modified time: 2021-12-17 23:36:30
import os

import yaml

# TODO: 这里需要将配置文件导入路径和 Pipe 类统一
MODULEPATH = os.path.dirname(os.path.abspath(__file__))
D_CONFIG = os.path.abspath(os.path.join(MODULEPATH, "../../config"))
_software = yaml.safe_load(open(f"{D_CONFIG}/software.yml", 'r').read())
_db = yaml.safe_load(open(f"{D_CONFIG}/database.yml", 'r').read())


def check_file(file_path):
    """
    Check whether a file exists

    :param file_path: The path of file
    :return: True or False
    """
    file_path = os.path.abspath(file_path)
    if os.path.isfile(file_path) and os.path.exists(file_path):
        return True
    else:
        return False


def check_dir(dir_path):
    """
    Check whether a dir exists

    :param dir_path: The path of dir
    :return: True of False
    """
    dir_path = os.path.abspath(dir_path)
    if os.path.isdir(dir_path) and os.path.exists(dir_path):
        return True
    else:
        return False


def single_run(script):
    """
    Run the script with shell

    :params script: The script path you want to exec
    """
    return f"{_software['shell']} {script}"


def multi_run(script, maxjob=4):
    """
    Use parallel to run the script

    :params script: The script path you want to exec
    :params maxjob:  The maxjob you want to exec at the same time
    """
    return f"{_software['parallel']} -j {maxjob} < {script}"


def parallel_run(script, maxjob=4):
    """
    Use parallel to run the script

    :params script: The script path you want to exec
    :params maxjob:  The maxjob you want to exec at the same time
    """
    return f"{_software['parallel']} -j {maxjob} < {script}"


def qsub_run(script, maxjob=4, queue="defaultApp", jobprefix="work", resource="1G"):
    """
    Run the multi line of script at cluster

    :params script: The script path you want to exec
    :params queue:  The queue you want to use, default all.q
    :params maxjob:  The maxjob you want to exec at the same time
    :params jobprefix:  The job prefix to use
    :params resource:  The memory use for one job
    """
    return f"{_software['perl']} {_software['qsub']} --queue {queue} -maxjob {maxjob} --jobprefix {jobprefix} --resource vf={resource} {script}"


def CLASS(*args):
    """
    class is a reserved word in Python

    :param args:
    :return:
    """
    return {"class": ' '.join(args)}
