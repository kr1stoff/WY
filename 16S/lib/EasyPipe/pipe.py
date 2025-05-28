#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: MingJia
# @Date:   2020-10-27 10:10:51
# @Last Modified by:   MingJia
# @Last Modified time: 2021-12-17 10:10:51
import subprocess
import logging
import os
import yaml

from .jobs import ComplexJob
from .jobs import Job
from .utils import check_dir
from .utils import check_file
from .utils import multi_run
from .utils import single_run

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

MODULEPATH = os.path.dirname(os.path.abspath(__file__))


class Pipe(ComplexJob):
    """A Class for building pipeline easy.

    1. Import the module, begain an object.
    2. Start the pipe with the method step_start.
    3. Use the method write_command to write command to the step
    4. End this step with method step_end
    You could repeat 2~4 to generate more step.

    ```
    import sys
    sys.path.append(os.path.join('/home/renchaobo/TestDir'))
    from EasyPipe import Pipe

    project = Pipe('GDD0000','./Pipe')
    qc = project.child("QC",f"{project.out}/QC",_type = "Job")
    qc.add_command("This is QC Step")
    align = project.child("align",f"{project.out}/Align")
    bw2 = align.child("BW2",f"{align.out}/BW2",_type = "Job")
    bw2.add_command("bw2_1")
    bw2.add_command("bw2_2")
    bw2.add_command("bw2_3")
    bw2stat = align.child("bw2stat",f"{align.out}/Align",_type = "Job")
    bw2stat.add_command("This is Align stat")
    project.finish()
    ```
    """

    def __init__(self,
                 name,
                 out,
                 bin_folder=None,
                 config_folder=None):
        """
        Init the pipe.
        You must offer the name and the project_dir.
        This will make the dirs need for the pipeline.

        :param name: The project name
        :param out: The out put dir
        :param bin_folder: The default bin dir
        :param config_folder: The default config dir contain
        """
        super(Pipe, self).__init__(name, out)

        # some private variable
        self.name = name
        self.out = os.path.abspath(out)
        os.makedirs(self.out, exist_ok=True)

        # 脚本目录
        self.script = os.path.join(self.out, 'SH')
        os.makedirs(self.script, exist_ok=True)
        # 日志
        self.log = os.path.join(self.script, "run.log")
        self.err = os.path.join(self.script, "run.err")
        self.workdirs = {}

        ## 存放可执行脚本的目录
        self.set_bin(bin_folder)
        ## 存放配置文件的目录
        self._d_config = os.path.abspath(config_folder) if config_folder else f"{MODULEPATH}/config"

        # 信息，软件，环境，数据库，bin
        # self.infos 可以用来存储一些分组信息，差异信息内容
        # 如果不提供包含database.yml及software.yml的config_folder参数，
        # 则从模块自带的配置文件种读取
        f_software = f"{self._d_config}/software.yml" if os.path.exists(
            f"{self._d_config}/software.yml") else f"{MODULEPATH}/config/software.yml"
        f_database = f"{self._d_config}/database.yml" if os.path.exists(
            f"{self._d_config}/database.yml") else f"{MODULEPATH}/config/database.yml"
        self.software = yaml.safe_load(open(f"{f_software}", 'r').read())
        self.db = yaml.safe_load(open(f"{f_database}", 'r').read())

    def __str__(self):
        res = f"""ProjectName : {self.name}
ProjectDir  : {self.out}
BinDir      : {self.BIN}
Softwares   : {self.software}
Databases   : {self.db}
"""
        return res

    def finish(self, run=0):
        """
        Output the scripts

        :param run: Whether run the pipe
        """
        script_all = os.path.join(self.script,
                                  f"{self.name}_All.sh")
        handle = open(script_all, 'w')
        handle.write(f"# {self.name}\nset -e\n")
        i = 0
        for step, element in self.children.items():
            handle.write(f'echo "Start {element.name} at $(date +%Y-%m-%d\ %H:%M:%S)"\n')
            i += 1
            step_script_name = os.path.join(self.script, f"Step{i}_{step}.sh")
            step_script_dir = os.path.join(self.script, f"Step{i}")

            if isinstance(element, Job):
                step_script_handle = open(step_script_name, 'w')
                sub_step_script = os.path.join(step_script_dir,
                                               f"{element.name}.sh")
                element.to_script(sub_step_script)
                if element.params["run_type"] == "single_run":
                    step_script_handle.write(single_run(sub_step_script) + "\n")
                elif element.params["run_type"] == "multi_run":
                    step_script_handle.write(multi_run(sub_step_script,
                                                       maxjob=element.params["maxjob"]) + "\n")
                else:
                    raise ValueError(
                        f"run_type must be single_run, multi_run")
                step_script_handle.close()
            elif isinstance(element, ComplexJob):
                element.to_script(step_script_name, step_script_dir)
            else:
                raise ValueError(f"element must Job or ComplexJob obj")
            handle.write(single_run(step_script_name) + "\n")
            handle.write(f'echo "End {element.name} at $(date +%Y-%m-%d\ %H:%M:%S)"\n')

        handle.close()
        if run:
            subprocess.run(f"bash {script_all}", shell=True)

    def set_bin(self, bin_path):
        """
        Set the bin dir path
        The bin is the dir where store the scripts used in the your pipeline

        :return: True or False
        """
        bin_path = os.path.abspath(bin_path)
        if check_dir(bin_path):
            self.bin = bin_path
            return True
        else:
            logging.warning(f"NO DIR {bin_path}")
            return False

    def add_software(self, key, value, check=False):
        """
        Add softwares the pipe use
        """
        if check:
            if check_file(value):
                self.softwares[key] = value
            else:
                raise Exception(f"Software path Error: {value}")
        else:
            self.softwares[key] = value

    def add_db(self, key, value, check=False):
        """
        Add dbs the pipe use
        """
        if check:
            if check_file(value):
                self.dbs[key] = value
            else:
                raise Exception(f"Database path Error: {value}")
        else:
            self.dbs[key] = value
