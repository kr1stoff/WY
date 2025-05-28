#!/usr/bin/env python
# @Author: zhuzi
# @Date:   2023/5/10
import os
import logging
import yaml
import click
import lib.common as common
import sys
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
PATH = os.path.dirname(os.path.abspath(__file__))

params = common.config()   # 读取配置文件
p_src = os.path.join(PATH, "lib/report/Tree_Typing/src")


class MultiRun():
    def __init__(self,out,upload,kind):
        self.kind = kind
        self.out = os.path.abspath(out)
        self.upload = os.path.abspath(upload)
        self.p_sh = os.path.join(self.out,'0.shell')
        os.makedirs(self.p_sh,exist_ok=True)


    def single_denovo(self,table,ref):
        path_sh_list = []
        _yaml = yaml.safe_load(open(table, 'r').read())
        self.id_list = list(_yaml['samples'].keys())
        ref = 'None' if ref == 'None' else os.path.abspath(ref)

        for id in self.id_list:
            cmd = []
            sh = os.path.join(self.p_sh,f'run_{id}.sh')
            path_sh_list.append(sh)   #记录每个样本的shell脚本路径
            cmd.append(f"time {params.python3} {PATH}/main.py -i {id} -p {self.out}/{id}/pipeline.yaml -s {self.out}/task_samples.yaml -r {ref}")
            common.cmd2shell(cmd,sh)

            # 进化树类型
            if self.kind=='bacteria':  self.tree_mode='core'
            elif self.kind=='fungi' :  self.tree_mode='core'
            elif self.kind=='virus' :  self.tree_mode='wgs'
        return path_sh_list


    def prepare_fa(self):
        sh = os.path.join(self.p_sh,'cp_fa.sh')
        self.fa_dir = os.path.join(self.out,'fa')
        os.makedirs(self.fa_dir,exist_ok=True)
        cmd = []
        for id in self.id_list:
            if self.kind != 'virus':
                cmd.append(f"cp -r {self.out}/{id}/{id}/2.ass/{id}_scaffolds.fasta {self.fa_dir}/{id}.fa")
            else:
                cmd.append(f"cp -r {self.out}/{id}/{id}/2.ass/{id}_consensus.fa {self.fa_dir}/{id}.fa")
        common.cmd2shell(cmd,sh)
        return sh


    # 进化树
    def tree(self,d_tree):
        sh = os.path.join(self.p_sh,'tree.sh')
        p_tree = f"{self.out}/Tree"
        p_yaml = f"{p_tree}/tree.yaml"
        os.makedirs(p_tree,exist_ok=True)

        fa_dic={}
        for file in os.listdir(self.fa_dir):
            id=file.replace('.fa','')
            fa_dic[id] = os.path.join(self.fa_dir,file)

        for dir in os.listdir(d_tree):
            fa_dic[dir] = os.path.join(d_tree,dir)

        tree_dic = {}
        tree_dic['samples'] = {}
        tree_dic['samples'] = fa_dic
        tree_dic['library'] = 'Tree'
        tree_dic['result_dir'] = p_tree
        tree_dic['kindom'] = self.kind

        #生成yaml
        with open(p_yaml, "w", encoding="utf-8", newline="") as f:
            f.write(yaml.dump(tree_dic, allow_unicode=True))

        cmd = []
        cmd.append(f"{params.python3} {params.p_tree}/main.py -p {self.tree_mode} -i {p_yaml}")   #run进化树脚本
        for id in self.id_list: #拷贝结果到各样本的结果目录
            cmd.append(f"mkdir -p {self.out}/{id}/{id}/Tree_Typing/Tree")
            cmd.append(f"cp -r {p_tree}/Tree/Upload/*  {self.out}/{id}/{id}/Tree_Typing/Tree")
        common.cmd2shell(cmd,sh)
        return sh


    # 分型
    def typing(self,typing):
        cmd = []
        sh = os.path.join(self.p_sh,'typing.sh')
        for id in self.id_list:
            p_typing = f"{self.out}/{id}/typing"
            p_res = f"{self.out}/{id}/{id}/Tree_Typing/Typing"
            cmd.append(f"""# typing
mkdir -p {p_typing}
mkdir -p {p_res}
{params.python3} {params.p_typing}/main.py --samples {id} -d {typing} --fastas {self.fa_dir}/{id}.fa -o {p_typing}
cp -r {p_typing}/Upload/{id} {p_res}
""")
        common.cmd2shell(cmd,sh)
        return sh


    # 报告
    def report(self):
        sh = os.path.join(self.p_sh ,'report.sh')
        cmd = []
        for id in self.id_list:
            cmd.append(f"""#report
cp -r {params.p_GDHR}/src {self.out}/{id}/{id}/Tree_Typing
{params.perl} {PATH}/lib/report/tree_typing_report.pl {id} {self.out}/{id}/{id}/Tree_Typing

cd {self.out}/{id}
zip -qr {id}.zip {id}
cd {self.out}/{id}/{id}
zip -qr Tree_Typing.zip Tree_Typing

cp -r {self.out}/{id}/{id} {self.upload}
cp -r {self.out}/{id}/{id}.zip {self.upload}
""")
        common.cmd2shell(cmd,sh)
        return sh




#### 参数 ================================================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-table', required=True,type=click.Path(),help="solar yaml文件")
@click.option('-kind', required=True,type=click.Choice(["bacteria","fungi","virus"]),help="病原体类型")
@click.option('-out', required=True,type=click.Path(),help="输出目录")
@click.option('-upload', required=False,type=click.Path(),default='None',help="结果上传目录")
@click.option('-ref', required=False,type=click.Path(),default='None',help="参考序列")
@click.option('-tree', required=False,default='blank',type=click.STRING,help="进化树目录")
@click.option('-typing', required=False,type=click.STRING,help="分型目录")

def main(table,kind,out,upload,tree,typing,ref):
    logfile=f"{out}/log"
    os.makedirs(logfile,exist_ok=True)

    project = MultiRun(out,upload,kind)
    
    # 单菌流程
    path_shell_list = project.single_denovo(table,ref)
    common.mul_pool(path_shell_list, logfile,params.jobs)   # 目前设置为 3 个

    sh2 = project.prepare_fa()
    common.prun((sh2, logfile))

    # 跑进化树
    if tree != "blank":
        tree = os.path.abspath(tree)
        sh3 = project.tree(tree)
        common.prun((sh3, logfile))

    # 跑分型
    if typing != "blank":   
        sh4 = project.typing(typing)
        common.prun((sh4, logfile))

    # 生成进化树/分型报告
    if tree != "blank" or  typing != "blank":
        sh5 = project.report()
        common.prun((sh5, logfile))


if __name__ == "__main__":
    main()
