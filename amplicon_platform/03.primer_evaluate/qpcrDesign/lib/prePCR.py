import logging
import subprocess
from pathlib import Path
import yaml
import re

# 日志格式
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class prePCR():
    def __init__(self, ref:str, bed:str, outdir:str, retry:int) -> None:
        self.ref = Path(ref).resolve()
        self.bed = Path(bed).resolve()
        self.outdir = Path(outdir).resolve()
        self.retry = retry
        self.config = self.parse_config_file()

    def parse_config_file(self):
        logging.info('获取配置文件信息.')
        yml = Path(__file__).parents[1].joinpath('conf/prepcr.yml')
        return yaml.safe_load(open(yml, 'r').read())

    def myrun(self, cmd:str):
        """运行命令行."""
        logging.info('运行命令: ' + re.sub(' +',' ',cmd))
        return subprocess.run(cmd, shell=True, executable='/bin/bash')

    def NGSPrimerPlex(self, bed:str, autoadjust=False):
        """运行NGS-PrimerPlex."""
        #复制bed
        self.myrun(f'ln -sf {bed} {self.outdir}/2tNGS.bed')
        #是否自动调整
        autoadj = '-autoadjust' if autoadjust else ''
        soft = self.config['soft']
        params = self.config['params']
        cmd = f"""
        {soft['python']} {soft['NGSPrimerPlex']} \
            -ref {self.ref} -region {self.outdir}/2tNGS.bed \
            -th {params['th']} -primernum1 {params['primernum1']} \
            -minampllen {params['minampllen']} -maxampllen {params['maxampllen']} \
            -optampllen {params['optampllen']} \
            -minprimerlen {params['minprimerlen']} -maxprimerlen {params['maxprimerlen']} \
            -optprimerlen {params['optprimerlen']} \
            -minprimermelt {params['minprimermelt']} -maxprimermelt {params['maxprimermelt']} \
            -optprimermelt {params['optprimermelt']} \
            -minprimergc {params['minprimergc']} -maxprimergc {params['maxprimergc']} \
            -optprimergc {params['optprimergc']} \
            -minprimerendgc {params['minprimerendgc']} -maxprimerendgc {params['maxprimerendgc']} \
            -maxprimerpolyn {params['maxprimerpolyn']}  -maxoverlap {params['maxoverlap']} \
            -returnvariantsnum {params['returnvariantsnum']}\
            -blast -skip {autoadj}
        """
        self.myrun(cmd)
    
    def xlsx2csv(self):
        """excel转csv格式."""
        soft = self.config['soft']
        cmd = f"""
        {soft['csvtk']} xlsx2csv -t {self.outdir}/2tNGS_primers_combination_1_info.xls \
            > {self.outdir}/2tNGS_primers_combination_1_info.txt
        """
        self.myrun(cmd)

    def is_plex_complete(self):
        """检查NGS-PrimerPlex运行成功,结果完整."""
        return Path(f'{self.outdir}/2tNGS_primers_combination_1_info.xls').resolve().exists()

    def get_names_from_log(self, log, bed, outbed):
        """NGS-PrimerPlex运行失败的log可以找到具体失败的bedname."""
        #拉丁名
        with open(bed) as f:
            line = next(f)
        latin = line.strip().split('\t')[3].split('-')[0]
        #discard bed name list
        dbns = []
        pat = re.compile(f'INFO.*    ({latin}-.*)_\d$')
        with open(log) as f:
            for line in f:
                res = re.findall(pat, line)
                if res:
                    dbns.append(res[0])
        #筛选bed文件
        with open(bed) as f, open(outbed, 'w') as g:
            for line in f:
                name = line.strip().split('\t')[3]
                if name not in dbns:
                    g.write(line)

    def execute(self):
        self.outdir.mkdir(parents=True, exist_ok=True)
        #第一次运行
        self.NGSPrimerPlex(self.bed)
        #检查是否成功,不成功加参数autoadjust
        if not self.is_plex_complete():
            self.NGSPrimerPlex(self.bed, autoadjust=True)
        #重试
        for i in range(self.retry - 1):
            #有结果就跳过
            if self.is_plex_complete(): 
                continue
            #第一次用self.bed, 第二次用上一次生成的
            if i:
                crntbed = f'{self.outdir}/{i-1}.bed'
            else:
                crntbed = self.bed #current bed
            self.get_names_from_log(f'{self.outdir}/2tNGS.log', crntbed, f'{self.outdir}/{i}.bed')
            self.NGSPrimerPlex(f'{self.outdir}/{i}.bed')
        #final xlsx2csv
        if self.is_plex_complete():
            self.xlsx2csv()
