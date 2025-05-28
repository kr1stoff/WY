from pathlib import Path
import re
from subprocess import run
import logging
import pandas as pd
from Bio.Seq import Seq
import sys
sys.path.append(Path(__file__).parent)
from .config import config



class DegenerateMSA():
    def __init__(self, primertable, allfna, outdir) -> None:
        self.params = config()
        self.primtbl = primertable
        self.allfna = allfna
        self.outdir = Path(outdir).resolve()

    def parse_merge_primer(self):
        """
        提取内左右引物/探针/外左右引物/参考登录号.
        return[0]:
        rec.array([(0, 'NC_001538.1', 'c60f46c32fae373f6c5801d8f67f7cd0', 'AGGAGCTGCTTACCCATGGAAT', 'AGGCTCGCAAAACATGTCTGT', 
                'CAGCCAAACCATGACCTCAGGAAGGA', 'GGAATGCAGCCAAACCATGA', 'TGGGGACAAGGCCAAGATT'),
                (1, 'NC_001538.1', 'cd99843686280449446e1d6cd580252e', 'CGGGACTGTAACACCTGCTCTT', 'GGGTTCCTTTGGCTTTTTGG', 
                'CCCCAACCAAAAGAAAAGGAGAGTGTCCA', 'CTATTGCCCCAGGAGGTGCTAAT', 'CCTCCTTTTATTAGTAGTTTTGGCACT'),
                ...
                dtype=[('index', '<i8'), ('Chrom', 'O'), ('ID', 'O'), ('Left_inSeq', 'O'), ('Right_inSeq', 'O'), ('Probe_Seq', 'O'), 
                ('Left_outSeq', 'O'), ('Right_outSeq', 'O')])
        """
        df_primtbl = pd.read_table(self.primtbl, sep='\t')
        #后续表头变化要对应修改!!
        df_chrom_seq = df_primtbl.loc[:, ['Chrom', 'ID', 'Left_Seq', 'Right_Seq', 'Probe_Seq', 'Left_Seq.1', 'Right_Seq.1']]
        df_chrom_seq.to_csv(self.outdir.joinpath('chrom_seq.tsv'), sep='\t', index=False)
        arr_chrom_seq = df_chrom_seq.to_records()
        return arr_chrom_seq, df_primtbl

    def cmds_to_sh(self, shlfl, cmds):
        """命令行写入Shell文件."""
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

    def get_amplicon_degenerate_seq(self, arr_chrom_seq):
        """根据外引物提取all.fna中扩增子序列, 检查扩增子中的简并和低频SNP"""
        allshl = []
        script_degenerate_check = Path(__file__).parents[1].joinpath('bin/degenerage_check.py')
        for rec in arr_chrom_seq:
            chrom, hsid, louts, routs = rec['Chrom'], rec['ID'], rec['Left_Seq.1'], rec['Right_Seq.1']
            cmd = f"""
            {self.params.seqkit} amplicon --bed -m 2 -j 32 -F {louts} -R {routs} {self.allfna} > {self.outdir}/{hsid}.bed
            awk -F '\\t' '{{printf ">%s\\n%s\\n",$1,$7}}' {self.outdir}/{hsid}.bed > {self.outdir}/{hsid}.fa
            {self.params.python} {script_degenerate_check} -s {self.outdir}/{hsid}.fa -n {chrom} -o {self.outdir}/{hsid}
            """
            cmd = re.sub(' +', ' ', cmd.strip())
            self.cmds_to_sh(f'{self.outdir}/{hsid}.sh', cmd)
            allshl.append(f'bash {self.outdir}/{hsid}.sh')
        shlfl = self.outdir.joinpath('seqkit_amplicon.sh')
        self.cmds_to_sh(shlfl, allshl)
        run(f'cat {shlfl} | {self.params.parallel} -j {self.params.prlnum}', shell=True, executable='/bin/bash')

    def degfunc(self, refbs, degbs, prim, orientation='F'):
        """
        根据引物在ref序列中的位置找在deg序列中找简并转化后的序列.
        orientation 引物方向
            F : forward primer
            R : reverse primer
            P : probe
        """
        if orientation == 'R':
            prim = str(Seq(prim).reverse_complement())
        res = re.search(prim, refbs)
        #检查有没有匹配到位置
        if res:
            span = res.span()
            if orientation == 'R':
                oseq = str(Seq(degbs[span[0]:span[1]]).reverse_complement())
            else:
                oseq = degbs[span[0]:span[1]]
        else:
            logging.error(f'没有匹配到参考上的位置.\n\tref:{refbs}\n\tquery:{prim}\n\torientation:{orientation}')
            sys.exit(1)
        return oseq

    def get_primer_degenerate_seq(self, dfp, arr_chrom_seq):
        """获取所有引物探针的简并转化序列."""
        for rec in arr_chrom_seq:
            chrom, hsid, louts, routs = rec['Chrom'], rec['ID'], rec['Left_Seq.1'], rec['Right_Seq.1']
            lins, rins, prbs = rec['Left_Seq'], rec['Right_Seq'], rec['Probe_Seq']
            with open(f'{self.outdir}/{hsid}/degenerate.txt') as f:
                refbs = next(f).strip()
                degbs = next(f).strip()
            dfp.loc[dfp['ID'] == hsid, 'Left_Seq'] = self.degfunc(refbs, degbs, lins)
            dfp.loc[dfp['ID'] == hsid, 'Right_Seq'] = self.degfunc(refbs, degbs, rins, orientation='R')
            dfp.loc[dfp['ID'] == hsid, 'Probe_Seq'] = self.degfunc(refbs, degbs, prbs)
            dfp.loc[dfp['ID'] == hsid, 'Left_Seq.1'] = self.degfunc(refbs, degbs, louts)
            dfp.loc[dfp['ID'] == hsid, 'Right_Seq.1'] = self.degfunc(refbs, degbs, routs, orientation='R')
        dfp.to_csv(f'{self.outdir}/{Path(self.primtbl).stem}.degenerate.xls', sep='\t', index=False)

    def degenerate(self):
        """整合简并流程."""
        self.outdir.mkdir(exist_ok=True, parents=True)
        arr_chrom_seq, df_primtbl = self.parse_merge_primer()
        self.get_amplicon_degenerate_seq(arr_chrom_seq)
        self.get_primer_degenerate_seq(df_primtbl, arr_chrom_seq)
