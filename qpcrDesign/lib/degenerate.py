from pathlib import Path
import re
from subprocess import run
import logging
import pandas as pd
import sys
from collections import namedtuple
import yaml


# 一些函数
def config():
    """ 读配置信息写到namedtuple中. e.g. params.samtools """
    config_yaml = Path(__file__).parents[1].joinpath('conf/config.yaml')
    dict_conf = yaml.safe_load(open(config_yaml))
    Params = namedtuple("Params", dict_conf.keys())
    params = Params(**dict_conf)
    return params


def cmds_to_sh(shlfl, cmds):
    """ 命令行写入Shell文件 """
    if isinstance(cmds, str):
        with open(shlfl, 'w') as g:
            g.write(cmds)
    elif isinstance(cmds, list):
        with open(shlfl, 'w') as g:
            for cmd in cmds:
                g.write(cmd + '\n')
    else:
        logging.error(f'不存在的命令行格式: {cmds}')
        sys.exit(1)


def get_reverse_complement(seq:str):
    """ 反向互补序列 """
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R', 'M': 'K', 
              'K': 'M', 'S': 'W', 'W': 'S', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N', 
              'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'r': 'y', 'y': 'r', 'm': 'k', 'k': 'm', 
              's': 'w', 'w': 's', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n'}
    seq_revcom = ''.join(reversed([comp_dict[s] for s in seq]))
    return seq_revcom


def find_dgnrt_seq_from_primers(refbs, degbs, prmbs, orientation='F'):
    """
    根据引物在ref序列中的位置找在deg序列中找简并转化后的序列.
    orientation 引物方向
        F : forward primer
        R : reverse primer
        P : probe
    """
    if orientation == 'R':
        prmbs = get_reverse_complement(prmbs)
    res = re.search(prmbs, refbs)
    #检查有没有匹配到位置
    if res:
        span = res.span()
        if orientation == 'R':
            oseq = get_reverse_complement(degbs[span[0]:span[1]])
        else:
            oseq = degbs[span[0]:span[1]]
    else:
        logging.error(f'没有匹配到参考上的位置.\n\tref:{refbs}\n\tquery:{prmbs}\n\torientation:{orientation}')
        sys.exit(1)
    return oseq


def degenerate_4base_depth_to_html(prmseq, degseq, dict_pos_4bs_depth, orientation='F'):
    """ 获取简并/小写碱基的4碱基深度信息, 直接写成HTML<a> """
    content = ''
    # 左引物
    if orientation == 'F':
        # 搜序列在简并序列中的起始位置
        seqstart = re.search(prmseq, degseq).span()[0]
        # 每个简并和小写的碱基都打印出来深度信息
        for idx in range(len(prmseq)):
            if prmseq[idx] not in 'ATGC':
                crrnt_pos = seqstart + idx
                content += f'<a title="{dict_pos_4bs_depth[crrnt_pos]}">{prmseq[idx]}</a>'
            else:
                content += prmseq[idx]
    # 右引物
    else:
        rcprmseq = get_reverse_complement(prmseq) # reverse complement primer sequence 
        seqstart = re.search(rcprmseq, degseq).span()[0]
        # 反向互补引物, 先把深度信息按顺序写到列表, 再反正取回来
        depths = []
        for idx in range(len(rcprmseq)):
            if rcprmseq[idx] not in 'ATGC':
                crrnt_pos = seqstart + idx
                A, T, G, C = dict_pos_4bs_depth[crrnt_pos].split(',')
                # 右引物的深度要转成互补碱基
                depths.append(','.join([T, A, C, G]))
        # 右引物顺序是反的
        depths.reverse()
        flag = 0 # 记录简并位置, 提取深度信息
        for idx in range(len(prmseq)):
            if prmseq[idx] not in 'ATGC':
                content += f'<a title="{depths[flag]}">{prmseq[idx]}</a>'
                flag += 1
            else:
                content += prmseq[idx]
    return content


################################################################################


class DegenerateMSA():
    def __init__(self, primertable, allfna, outdir) -> None:
        self.params = config()
        self.prmtbl = primertable
        self.allfna = allfna
        self.outdir = Path(outdir).resolve()

    def extract_info_from_merge_primer_table(self):
        """
        提取内左右引物/探针/外左右引物/参考登录号, 返回原始核心引物信息数组, 返回原始全部引物信息表
        return[0]:
        rec.array([(0, 'NC_001538.1', 'c60f46c32fae373f6c5801d8f67f7cd0', 'AGGAGCTGCTTACCCATGGAAT', 'AGGCTCGCAAAACATGTCTGT', 
                'CAGCCAAACCATGACCTCAGGAAGGA', 'GGAATGCAGCCAAACCATGA', 'TGGGGACAAGGCCAAGATT'),
                ...
                dtype=[('index', '<i8'), ('Chrom', 'O'), ('ID', 'O'), ('Left_inSeq', 'O'), ('Right_inSeq', 'O'), ('Probe_Seq', 'O'), 
                ('Left_outSeq', 'O'), ('Right_outSeq', 'O')])
        """
        df_prm_tbl_org = pd.read_table(self.prmtbl, sep='\t')
        # 后续表头变化要对应修改!!
        df_prm_info_org = df_prm_tbl_org.loc[:, ['Chrom', 'ID', 'Left_Seq', 'Right_Seq', 'Probe_Seq', 'Left_Seq.1', 'Right_Seq.1']]
        # df_prm_info_org.to_csv(self.outdir.joinpath('prm_info.tsv'), sep='\t', index=False)
        arr_prm_info_org = df_prm_info_org.to_records()
        return arr_prm_info_org, df_prm_tbl_org

    def batch_analyze_amplicon_degenerate(self):
        """ 根据外引物提取all.fna中扩增子序列, 分析扩增子中的简并和低频SNP """
        allshl = []
        script_degenerate_check = Path(__file__).parents[1].joinpath('bin/degenerage_analysis.py')
        for rec in self.arr_prm_info_org:
            chrom, hsid, louts, routs = rec['Chrom'], rec['ID'], rec['Left_Seq.1'], rec['Right_Seq.1']
            cmd = f"""
{self.params.seqkit} amplicon --bed -m 2 -j 32 -F {louts} -R {routs} {self.allfna} > {self.outdir}/{hsid}.bed
awk -F '\\t' '{{printf ">%s\\n%s\\n",$1,$7}}' {self.outdir}/{hsid}.bed > {self.outdir}/{hsid}.fa
{self.params.python} {script_degenerate_check} -s {self.outdir}/{hsid}.fa -n {chrom} -o {self.outdir}/{hsid}
            """
            cmd = re.sub(' +', ' ', cmd.strip())
            cmds_to_sh(f'{self.outdir}/{hsid}.sh', cmd)
            allshl.append(f'bash {self.outdir}/{hsid}.sh')
        shlfl = self.outdir.joinpath('seqkit_amplicon.sh')
        cmds_to_sh(shlfl, allshl)
        run(f'cat {shlfl} | {self.params.parallel} -j {self.params.prlnum}', shell=True, executable='/bin/bash')

    def convert_to_degenerate_primer_table(self):
        """ 获取所有引物探针的简并转化序列, 返回简并转换后的核心引物信息表 """
        dfp = self.df_prm_tbl_org
        for rec in self.arr_prm_info_org:
            chrom, hsid, louts, routs = rec['Chrom'], rec['ID'], rec['Left_Seq.1'], rec['Right_Seq.1']
            lins, rins, prbs = rec['Left_Seq'], rec['Right_Seq'], rec['Probe_Seq']
            with open(f'{self.outdir}/{hsid}/degenerate.txt') as f:
                refbs = next(f).strip()
                degbs = next(f).strip()
            dfp.loc[dfp['ID'] == hsid, 'Left_Seq'] = find_dgnrt_seq_from_primers(refbs, degbs, lins)
            dfp.loc[dfp['ID'] == hsid, 'Right_Seq'] = find_dgnrt_seq_from_primers(refbs, degbs, rins, orientation='R')
            dfp.loc[dfp['ID'] == hsid, 'Probe_Seq'] = find_dgnrt_seq_from_primers(refbs, degbs, prbs)
            dfp.loc[dfp['ID'] == hsid, 'Left_Seq.1'] = find_dgnrt_seq_from_primers(refbs, degbs, louts)
            dfp.loc[dfp['ID'] == hsid, 'Right_Seq.1'] = find_dgnrt_seq_from_primers(refbs, degbs, routs, orientation='R')
        dfp.to_csv(f'{self.outdir}/{Path(self.prmtbl).stem}.degenerate.xls', sep='\t', index=False)
        df_prm_tbl_deg = dfp.loc[:, ['Chrom', 'ID', 'Left_Seq', 'Right_Seq', 'Probe_Seq', 'Left_Seq.1', 'Right_Seq.1']]
        arr_prm_info_deg = df_prm_tbl_deg.to_records()
        return arr_prm_info_deg

    def analyze_degenerate_base_component(self):
        """ 分析简并碱基组成, 生成HTML网页 """
        html_tmplt = Path(__file__).parents[1].joinpath('etc/degenerate_template.html')
        html_dgnrt = f'{self.outdir}/{Path(self.prmtbl).stem}.degenerate.html'
        with open(html_dgnrt, 'w') as g, open(html_tmplt) as f:
            # HTML 模板
            cntt_tmplt = f.read()
            tbl_rows = ''
            for rec in self.arr_prm_info_deg:
                chrom, hsid, louts, routs = rec['Chrom'], rec['ID'], rec['Left_Seq.1'], rec['Right_Seq.1']
                lins, rins, prbs = rec['Left_Seq'], rec['Right_Seq'], rec['Probe_Seq']
                # 当前扩增子的4碱基深度字典. {0: '0,0,0,149', 1: '0,149,0,0', 2: '0,149,0,0', 3: '0,149,0,0', ...}
                dict_pos_4bs_depth = {}
                with open(f'{self.outdir}/{hsid}/result.depth.csv') as f:
                    next(f)
                    for line in f:
                        lns = line.strip().split(',')
                        dict_pos_4bs_depth.setdefault(int(lns[0]), ','.join(lns[1:]))
                # 两行简并文件，第一行参考，第二行简并转换后的
                with open(f'{self.outdir}/{hsid}/degenerate.txt') as f:
                    next(f)
                    degseq = next(f).strip()
                tbl_rows += (f"""
        <tr>
          <td>{chrom}</td>
          <td>{hsid}</td>
          <td>{degenerate_4base_depth_to_html(louts, degseq, dict_pos_4bs_depth)}</td>
          <td>{degenerate_4base_depth_to_html(routs, degseq, dict_pos_4bs_depth, orientation='R')}</td>
          <td>{degenerate_4base_depth_to_html(lins, degseq, dict_pos_4bs_depth)}</td>
          <td>{degenerate_4base_depth_to_html(rins, degseq, dict_pos_4bs_depth, orientation='R')}</td>
          <td>{degenerate_4base_depth_to_html(prbs, degseq, dict_pos_4bs_depth)}</td>
        </tr>
                        """)
            cntt_dgnrt = cntt_tmplt.replace('{% block tbody %}', tbl_rows)
            g.write(cntt_dgnrt)

    def degenerate(self):
        """整合简并流程."""
        self.outdir.mkdir(exist_ok=True, parents=True)
        self.arr_prm_info_org, self.df_prm_tbl_org = self.extract_info_from_merge_primer_table()
        self.batch_analyze_amplicon_degenerate()
        self.arr_prm_info_deg = self.convert_to_degenerate_primer_table()
        self.analyze_degenerate_base_component()
