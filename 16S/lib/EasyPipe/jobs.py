#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: MingJia
# @Date:   2020-10-27 10:10:51
# @Last Modified by:   MingJia
# @Last Modified time: 2021-12-17 10:10:51
import logging
import os

from .utils import multi_run
from .utils import single_run

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class Job(object):
    """
    """

    def __init__(self, name, out, run_type="single_run", maxjob=1):
        """
        Init the class

        :param name: The job name
        :param out: The out put dir
        :param run_type: The run type for the job(single_run|multi_run),default is single_run
        :param maxjob: If the run_type is multi_run, set the parallel run job number
        """
        self.name = name
        self.workdir = self.set_workdir(out)
        self.commands = []
        # 存储输出文件/目录
        self.output = {}
        self.params = {"run_type": run_type,
                       "maxjob": maxjob}
        if run_type == "single_run":
            self.set_maxjob(1)

    def add_command(self, command=None):
        """
        Add command content to the command list

        :param command: The command contetn you want to add to job
        """
        self.commands.append(command)

    @property
    def command(self):
        """
        Get the command content
        """
        return '\n'.join(self.commands)

    def set_workdir(self, workdir):
        """
        Set the out put result dir for the job
        """
        os.makedirs(workdir, exist_ok=True)
        self.workdir = workdir

    def set_run_type(self, run_type):
        """
        Set the run type for the complex job(single_run | multi_run)
        """
        if run_type not in set(["single_run", "multi_run"]):
            raise ValueError(
                f"run_type must be single_run, multi_run")
        else:
            self.params["run_type"] = run_type

    def set_maxjob(self, maxjob):
        self.params["maxjob"] = maxjob

    def to_script(self, script):
        """
        Write the command to a script file
        """
        dir_name = os.path.dirname(script)
        os.makedirs(dir_name, exist_ok=True)
        if self.params["maxjob"] == 1:
            self.commands = [f"set -e\n# {self.name}"] + self.commands
        with open(script, 'w') as OUT:
            print('\n'.join(self.commands), sep="\n", file=OUT)


class ComplexJob(object):
    """
    The ComplexJob class represent the complex jobs
    """

    def __init__(self, name, out, run_type="single_run", maxjob=1):
        """
        Init the class

        :param name: The job name
        :param out: The out put dir
        :param run_type: The run type for the job(single_run|multi_run),default is single_run
        :param maxjob: The max job number at the same time

        """
        self.name = name
        self.workdir = self.set_workdir(out)
        self.command = []
        self.children = {}
        # 存储输出文件/目录
        self.output = {}
        self.params = {"run_type": run_type,
                       "maxjob": maxjob}
        if run_type == "single_run":
            self.set_maxjob(1)

    def set_run_type(self, run_type):
        """
        Set the run type for the complex job(single_run | multi_run)
        """
        if run_type not in set(["single_run", "multi_run"]):
            raise ValueError(
                f"run_type must be single_run, multi_run")
        else:
            self.params["run_type"] = run_type

    def set_maxjob(self, maxjob):
        self.params["maxjob"] = maxjob

    def child(self, name, out, run_type="single_run", _type="ComplexJob"):
        """
        Generate a child

        :param name: The name for the child job/complexjob
        :param out: The out put dir for the child
        :param run_type: The run type for the job(single_run | multi_run)
        :param _type: The job type(Job | ComplexJob), default is ComplexJob
        """
        if _type == "Job":
            res = Job(name, out, run_type=run_type)
        elif _type == "ComplexJob":
            res = ComplexJob(name, out, run_type=run_type)
        self._add(res)
        return res

    def _add(self, element):
        """
        Add sub element to the object, it may be a Job object or a ComplexJob
        object

        :param element: Job object or ComplexJob object
        """
        assert (isinstance(element, Job) or isinstance(element,
                                                       ComplexJob)), f"Wrong Data Type for {element}"
        self.children[element.name] = element

    def set_workdir(self, workdir):
        """
        Set the work dir for the job

        :params workdir: The ComplexJob's workdir
        """
        os.makedirs(workdir, exist_ok=True)
        self.workdir = workdir

    def to_script(self, script, sub_script_dir):
        """
        Out put the run script

        :params script: The out put script to run
        :params sub_script_dir: The dir to put the sub scripts
        """

        os.makedirs(sub_script_dir, exist_ok=True)
        with open(script, 'w') as OUT:
            for name, element in self.children.items():
                if isinstance(element, Job):
                    script_name = os.path.join(sub_script_dir, f"{name}.sh")
                    element.to_script(script_name)
                    if element.params["run_type"] == "single_run":
                        OUT.write(single_run(script_name) + "\n")
                    elif element.params["run_type"] == "multi_run":
                        OUT.write(multi_run(script_name, maxjob=element.params["maxjob"]) + "\n")
                    else:
                        raise ValueError(
                            f"run_type must be single_run, multi_run or shell_run")
                elif isinstance(element, ComplexJob):
                    script_name = os.path.join(sub_script_dir, f"{name}.sh")
                    element.to_script(script_name, sub_script_dir)
                else:
                    raise ValueError(f"element must Job or ComplexJob obj")

    def _set_multi_run_type(self, run_type):
        """
        Set the multi run type

        :param run_type: The multi run type name (parallel_run|qsub_run)
        """
        if run_type not in {"parallel_run", "qsub_run"}:
            raise ValueError(f"Unknown multi_run type: {run_type}")
        else:
            self.multi_run_type = run_type

    def auto_multi_run_type(self):
        """
        Set the multi run type depend on the host name
        :return:
        """
        import platform
        hostname = platform.node()
        if hostname == "login":
            # GDIM
            self._set_multi_run_type("qsub_run")
        elif hostname == "localhost.localdomain":
            # home
            self._set_multi_run_type("parallel_run")
        elif hostname == "MZ72":
            # vision
            self._set_multi_run_type("parallel_run")
        else:
            raise ValueError(logging.error(f"Unknown HOST: {hostname}"))
