import os
import sys
import re
import logging
import subprocess
from multiprocessing import Pool
import time
from functools import wraps
import argparse
import glob
import zipfile
import shutil


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler)

################################################################################
#装饰器##########################################################################

def time_wrapper(func):
    """程序运行时间统计装饰器"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        r = func(*args, **kwargs)
        runtime = format_runtime(time.perf_counter() - start)
        logging.info(f"{func.__module__}.{func.__name__} : {runtime}")
        return r
    return wrapper

################################################################################
#运行类##########################################################################

class MyLoggingInfo():
    """让logging.info带上流程顺序自增标号"""
    def __init__(self):
        self.count = 0
    
    def info(self, loginfo):
        self.count += 1
        logging.info(f"{self.count}.{loginfo}")


class MyRule():
    """
    模仿snakemake工作流形式.
    [software]  软件
    [infile]    输入文件
    [log]       该软件的标准输出/标准错误的记录文件
    [ptn]       命令行通配字符串, e.g. "{software} -t {thread} -i {input}"
    [outfile]   (可选)输出文件     
    [params]    (可选)软件参数
    [thread]    (可选)软件线程数
    [runbool]   是否执行, 默认执行(True)
    """
    def __init__(self, software, infile, log, ptn, thread=1, params=None, outfile=None, runbool=True):
        self.cml = ptn.format(
            software=software,
            infile=infile,
            outfile=outfile,
            params=params,
            thread=thread
        )
        self.infile = infile
        self.outfile = outfile
        self.log = log
        # logging.info("shell: " + self.cml)
        myrun((self.cml, self.log))

class MyRuleS():
    """
    MyRule类的多样本并行模式, 参数与MyRule几乎一致, 但输入/输出/记录由字符串换为列表.
    注意infiles,logs,output需要列表长度一致并是相同样本!!!
    [infiles]   输入文件列表
    [logs]      记录文件列表
    [outfile]   (可选)输出文件
    ...
    其他参数同MyRule类
    """
    def __init__(self, software, infiles, logs, ptn, parallel=1, thread=1, params=None, 
    outfiles=None):
        self.software   = software
        self.infiles    = infiles
        self.logs       = logs
        self.ptn        = ptn
        self.params     = params
        self.outfiles   = outfiles
        self.thread     = thread
        self.parallel   = parallel
        self.generate_commands()
        self.run_multi()

    def generate_commands(self):
        """生成命令和记录文件的元组列表"""
        _infiles, _outputs, _logs = format_args_myrules(self.infiles, self.outfiles, self.logs)
        lists_command_log = list()
        # zitems: [0]infile [1]outfile [2]log
        for zitems in zip(_infiles, _outputs, _logs):
            cml = self.ptn.format(
                software=self.software,
                infile=zitems[0],
                outfile=zitems[1],
                params=self.params,
                thread=self.thread
            )
            # logging.info("shell: " + cml)
            lists_command_log.append((cml, zitems[2]))
        self.lists_command_log = lists_command_log
        
    def run_multi(self):
        with Pool(int(self.parallel)) as ph:
            ph.map(myrun, self.lists_command_log)


def format_args_myrules(*args):
    """MyRuleS, 使用该函数.功能: 把None格式化成[None,None,None...]列表,方便zip时候对齐"""
    arg_lengths = [len(arg) for arg in args if arg != None]
    args_out = ([None]*min(arg_lengths) if arg == None else arg for arg in args)
    return args_out


@time_wrapper
def myrun(command_log:tuple):
    """
    MyRule使用的运行函数,接收命令和记录文件的元组列表,兼容单样本和并行.
    加入环境获取,并把软件的BIN目录放在PATH第一位.
    [command_log] '(cml, log)' 命令行和记录文件组成的元组
    """
    _cml, _log = command_log
    # 环境
    dict_env = dict(os.environ)
    bin_current = os.path.dirname(_cml.split(" ")[0])
    dict_env["PATH"] = f"{bin_current}:{dict_env['PATH']}"
    # 运行
    res = subprocess.run(_cml, shell=True, encoding="utf-8",
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=dict_env)
    with open(_log, "at", encoding="utf-8", newline="") as gh:
        gh.write("[STDOUT]\n" + res.stdout)
        gh.write("\n[STDERR]\n" + res.stderr)

################################################################################
#综合函数#######################################################################

def get_args():
    """溯源进化树 - 命令行解析器"""
    parser = argparse.ArgumentParser(
        description='溯源进化树主程序, 包括全基因组多序列比对(MSA)进化树,SNP进化树和核心基因进化树.'
    )
    parser.add_argument('-p', '--pipe', choices=['wgs','snp','core'], required=True, 
                        help='选择要跑的子流程')
    parser.add_argument('-i', '--inyaml', help='样本信息YAML文件.', required=True)
    parser.add_argument('-n', '--dryrun', action='store_true', help='不运行,只生成脚本和目录. (default:False)')
    args = parser.parse_args()
    return args

def format_runtime(second):
    """秒转时分秒"""
    m, s = divmod(second, 60)
    h, m = divmod(m, 60)
    if h:
        return f"{h}h{m}m{s:.2f}s"
    elif m:
        return f"{m}m{s:.2f}s"
    else:
        return f"{s:.2f}s"


@time_wrapper
def merge_wgs_fasta(mfasta, dict_sample):
    """
    合并全基因组多序列比对FASTA
    [220801] customer tree 库文件名一样,里面的表头重复,什么鬼???格式化一下
    mfasta        : 输出合成的 FASTA
    dict_sample   : 样本字典,key是样本名, value是样本绝对路径
    """
    dir_mergefa = os.path.dirname(mfasta)
    cml = "cat "
    for name in dict_sample:
        fa = dict_sample[name]
        sub_fa = f"{dir_mergefa}/{os.path.basename(fa)}"
        sub_name = re.sub(".fasta|.fa|.fna", "", name)
        cml_sub = f"sed 's/>.*/>{sub_name}/g' {fa} > {sub_fa}"
        os.system(cml_sub)
        cml += f"{sub_fa} "
    cml += f"> {mfasta}"
    logging.debug(f"shell: {cml}")
    subprocess.run(cml, shell=True)


def format_fasta(name, infile, outfile):
    """
    格式化表头名
    [name]      : 将要格式化的样本名  
    [input]     : 输入 fasta
    [outfile]   : 输出 fasta
    """
    head_list   = list()
    headnum     = 1
    with open(infile, "rt") as fh, open(outfile, "wt", encoding="utf-8", newline="") as gh:
        for line in fh:
            if line.startswith(">"):
                newhead = f"{name}"
                while newhead in head_list:
                    headnum += 1
                    newhead = f"{name}_{headnum}"
                head_list.append(newhead)
                gh.write(f">{newhead}\n")
            else:
                gh.write(line)


def mymakedirs(dirs:list, outdir):
    """根据列表创建子目录"""
    for dr in dirs:
        os.makedirs(os.path.join(outdir, dr), exist_ok=True)


def svg2png(path_magick, svg):
    """
    svg图片转png
    """
    MyRule(
        software=path_magick,
        infile=svg,
        outfile=svg.replace("svg", "png"),
        params="-density 300 -colorspace RGB",
        log="/dev/null",
        ptn="{software} {params} {infile} {outfile}"
    )


@time_wrapper
def magick_dir(dir_fig, path_magick, suffix_from, suffix_to):
    """
    整个目录下的图片转换格式
    [dir_fig]       图片目录
    [path_magick]   magick软件路径
    [suffix_from]   模板图片后缀
    [suffix_to]     需要生成图片后缀
    """
    svgs = glob.glob(f"{dir_fig}/*.{suffix_from}")
    if suffix_from == "svg" and suffix_to == "png":
        for svg in svgs:
            svg2png(path_magick, svg)
    else:
        print("mxf: 添加判断.")


@time_wrapper
def zip_dir(dirpath, zip_dst):
    """压缩目录,和'linux zip -r'一样,不用切目录"""
    zip = zipfile.ZipFile(zip_dst, "w", zipfile.ZIP_DEFLATED)
    for root, dirs, files in os.walk(dirpath):
        arcpath = root.replace(os.path.dirname(dirpath), "")
        for file in files:
            zip.write(os.path.join(root, file), arcname=os.path.join(arcpath, file))
    zip.close()


def link_exist(src, dst):
    """
    只复制存在的文件, 并且没创建过
    [src] sourse源文件
    [dst] destination目的文件
    """
    if os.path.isfile(src) and not os.path.exists(dst):
        os.symlink(src, dst)


def link2upload(src, dir_upload):
    """源文件软连接到Upload目录"""
    link_exist(src, f"{dir_upload}/{os.path.basename(src)}")


def delete_wildcard_path(wildcard:str):
    """接收linux命令行一样的通配符删除文件"""
    list_res = glob.glob(wildcard)
    for path in list_res:
        if os.path.isdir(path):
            shutil.rmtree(path)
        elif os.path.isfile(path) or os.path.islink(path):
            os.remove(path)
        else:
            logging.warning(f"{path}不是文件夹也不是文件?")
