#!/usr/bin/env python
import os
import logging
from multiprocessing import Pool
from subprocess import Popen,PIPE,TimeoutExpired,run
import time
import signal
import yaml
from pathlib import Path
from collections import namedtuple
import sys

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)

def cmd2shell(cmds, sh):
    """
    [功能]：
        命令list写到shell脚本里
        cmds:  命令列表
        shell: 输出脚本
    """
    with open(sh, "w", encoding="utf-8", newline="") as w:
        w.write(f"#!/usr/bin/bash \n")
        w.write(f"echo -e [`date +%Y-%m-%d\" \"%T`] {os.path.basename(sh)} Start \n")
        for cmd in cmds:
            w.write(cmd + "\n")
        w.write(f"echo -e [`date +%Y-%m-%d\" \"%T`] {os.path.basename(sh)} End \n\n")
    os.system(f"chmod 755 {sh}")



def pe_fq_link(id,fq1,fq2,out):
    """
    [功能]：
        创建双端fq软连接
    """
    s_fix=str(os.path.basename(fq1))
    if s_fix.endswith(".gz"):
        l_fq1, l_fq2=f"{out}/{id}_1.fastq.gz" ,f"{out}/{id}_2.fastq.gz"
        creat_link(fq1,l_fq1)
        creat_link(fq2,l_fq2)
    elif s_fix.endswith((".fq",".fastq")):
        l_fq1, l_fq2=f"{out}/{id}_1.fastq" ,f"{out}/{id}_2.fastq"
        creat_link(fq1,l_fq1)
        creat_link(fq2,l_fq2)
    else:
        logging.debug(f"Fastq's suffix can't read ----({s_fix}) ")
        l_fq1 ,l_fq2=f"{out}/{id}_1.fastq.gz" ,f"{out}/{id}_1.fastq.gz"
        creat_link(fq1,l_fq1)
        creat_link(fq2,l_fq2)
    return l_fq1 ,l_fq2



def se_fq_link(id,fq,out):
    """
    [功能]：
        创建单端fq软连接
    """
    s_fix=str(os.path.basename(fq))
    if s_fix.endswith(".gz"):
        l_fq=f"{out}/{id}.fastq.gz" 
        creat_link(fq,l_fq)
    elif s_fix.endswith((".fq",".fastq")):
        l_fq=f"{out}/{id}.fastq" 
        creat_link(fq,l_fq)
    else:
        logging.debug(f"Fastq's suffix can't read ----({s_fix}) ")
        l_fq =f"{out}/{id}.fastq.gz" 
        creat_link(fq,l_fq)
    return l_fq



def creat_link(p_file, p_link):
    """
    [功能]：
        功能函数，创建软连接
    """
    if os.path.isfile(p_file) and not os.path.exists(p_link):
        os.symlink(p_file, p_link)
    else:
        logging.debug(f"不存在文件<{p_file}> 或 链接已存在<{p_link}>")


def kill_command(p):
    """ 杀进程 """
    os.killpg(os.getpgid(p.pid),signal.SIGTERM)


def prun(sh_tuple):
    """
    [功能]：
        单任务运行/进程池 函数
        sh_tuple = (shell命令文件 ，日志文件,超时时间)
    """
    shfile = sh_tuple[0]
    logdir = sh_tuple[1]
    if len(sh_tuple) == 3:
        timeout = sh_tuple[2]
    else:
        timeout = 5400

    shname=os.path.basename(shfile)
    O = open(f"{logdir}/{shname}.o",'w',encoding='utf-8')
    E = open(f"{logdir}/{shname}.e",'w',encoding='utf-8')

    try:
        process = Popen(shfile, stdout=PIPE, stderr=PIPE, shell=True,encoding='utf-8')
        process_id = process.pid    # 获取进程号
        stdout,stderr = process.communicate()
        O.write(str(stdout))
        E.write(str(stderr))

        start_time = time.time()
        while process.poll() is None:
            time.sleep(0.1)
            elapsed_time = time.time() - start_time # 时间差
            if elapsed_time > timeout:
                process.terminate()         # 如果超时，杀掉进程
                time.sleep(1)               # 等待一段时间，确保进程被终止
                E.write(f"### 任务超时阈值 {timeout} 秒，超时杀掉 {shfile} 引发的进程 {process_id}\n")
                
                # 如果进程仍然没有终止，强制终止
                if process.poll() is None:process.kill()

    except Exception as e:
        E.write(str(e))

    O.close()
    E.close()



def mul_pool(plist, logdir, num=None,timeout=5400):
    """
    [功能]：
        多任务并发 函数
        plist：  shell命令列表
        logdir： 日志文件
        num：    并发数量
        timeout：进程运行超时时间阈值
    """
    pool= Pool(num) if num else Pool(len(plist))
    argv_list=[(tmp,logdir,timeout) for tmp in plist]
    pool.map(prun, argv_list)



def get_cfg():
    """
    [功能]：
        获取系统的cpu / 线程数
    """
    # CPU 核数
    cmd = "cat /proc/cpuinfo| grep 'cpu cores'| uniq |cut -d' ' -f3"
    cpu = int(os.popen(cmd).read())

    # 线程数
    cmd = "cat /proc/cpuinfo| grep 'processor'| wc -l"
    threads = int(os.popen(cmd).read())
    return cpu,threads


def config():
    """读配置信息写到namedtuple中"""
    config_yaml = Path(__file__).parents[1].joinpath('cfg/config.yaml')
    dict_conf = yaml.safe_load(open(config_yaml))
    Params = namedtuple("Params", dict_conf.keys())
    params = Params(**dict_conf)
    return params


def dic2yaml(dic: dict, save_path: str):
    """dict保存为yaml"""
    with open(save_path, 'w' ,encoding='utf-8') as file:
        file.write(yaml.dump(dic, allow_unicode=True))
    logging.info(f"[Output]: dic===>yaml {save_path}")



def yaml2dic(yaml_path: str):
    """读取yaml"""
    with open(yaml_path,'r',encoding='utf-8') as file:
        dic = yaml.load(file.read(), Loader=yaml.FullLoader)
        return dic