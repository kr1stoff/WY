#!/usr/bin/env python3
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2022/9/22 16:40
import os
import re
import multiprocessing
import click
import logging
import sys

logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)



def run(cmd):
    print(cmd)
    os.system(cmd)

class island:
    def __init__(self,res):
        self.res = os.path.abspath(res)
        self.cmd_list = []

    def split_gbk(self,gbkFile):
        """
        Islandpath don't support set of gbk.
        So need to split a individual gbk file one by one
        """
        self.split_gbk_file = os.path.join(self.res,'Split_Gbk')
        os.makedirs(self.split_gbk_file,exist_ok=True)

        with open(gbkFile ,'r',encoding='utf-8') as f:
            for line in f:
                if not line: continue
                if re.search('^LOCUS' ,line):
                    li = re.split('\s+',line)
                    id = li[1]
                    tail = li[-3] +' '+ li[-2]
                    tmp = re.findall(r"(.*)_length_(\d+)_" ,id)[0]
                    short_id = tmp[0]
                    length = tmp[1]
                    OUT = open(f"{self.split_gbk_file}/{short_id}.gbk" ,'w' ,encoding='utf-8')
                    OUT.write(f"LOCUS       {short_id} {length} bp   {tail}  ")
                else:
                    OUT.write(line)
        OUT.close()


    def islandpath_sh(self,soft):
        """Write run island shell"""
        self.d_res = os.path.abspath(os.path.join(self.res,'result'))
        os.makedirs(self.d_res,exist_ok=True)
        os.chdir(self.res)

        sh = os.path.join(self.res,'work.sh')
        with open(sh ,'w' ,encoding='utf-8') as w:
            for file in os.listdir(self.split_gbk_file):
                if not file.endswith('.gbk'): continue
                id = file.replace('.gbk','')

                gbk_file = os.path.join(self.split_gbk_file,file)
                p_res = os.path.join(self.d_res,id)

                cmd = f"{soft} {gbk_file} {p_res}" 
                w.write(f"{cmd}\n")
                self.cmd_list.append(cmd)
        return self.cmd_list



    def parse_result(self):
        """grep need content"""
        cmd = f"""
        cd {self.d_res}
        grep 'gi' * |tr ':' '\\t' |sed 's/ID=unknown_//g' |awk '{{print $1\"\\t\"$10\"\\t\"$5\"\\t\"$6}}' |sort -Vk 1 >{self.res}/island.txt
        sed -i '1i Contig\\tGI_num\\tStart\\tEnd' {self.res}/island.txt
        """
        os.system(cmd)


#参数
p_soft = "/home/earthtest/miniconda3/envs/denovo/bin/islandpath"
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--gbk',required=True,type=click.Path(),help="Prokka gbk file")
@click.option('-o', '--out',required=True,type=click.Path(),help="Output dir")
@click.option('-p', '--path',type=click.Path(),show_default=True,default=p_soft,help="islandpath software path")

def main(gbk,out,path):
    """
    1.This pipeline need in denovo (conda'env)\n
    2.This pipeline need contig_headers such as [NODE_1_length_3451_cov_22]
    """

    project = island(out)
    project.split_gbk(gbk)
    cmd_list = project.islandpath_sh(path)

    pool = multiprocessing.Pool(processes=6)
    pool.map(run, cmd_list)

    project.parse_result()




#调用全局函数
if __name__ == "__main__":
    main()