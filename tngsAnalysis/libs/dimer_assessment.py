import logging
import re
from pathlib import Path
from collections import namedtuple, defaultdict
from subprocess import run
import pandas as pd
import matplotlib.pyplot as plt
from modules.config import config


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class DimerAssess():
    def __init__(self, fastq, primer_fa, outdir, prefix):
        """单端测序去引物"""
        self.fastq = fastq
        self.primer_fa = primer_fa
        self.outdir = outdir
        self.prefix = prefix
        #Others
        self.params = config()
        Path(self.outdir).joinpath('tmp').mkdir(parents=True, exist_ok=True)

    def seqtk_blastn(self):
        """
        1. seqtk fastq转fasta; 
        2. blast引物序列建库, blast比对reads和引物
        """
        primer_fa_name = Path(self.primer_fa).name
        cml = f"""
{self.params.seqtk} seq -A {self.fastq} > {self.outdir}/{self.prefix}.fasta
cp {self.primer_fa} {self.outdir}/tmp/{primer_fa_name}
{self.params.makeblastdb} -in {self.outdir}/tmp/{primer_fa_name} -dbtype nucl -out {self.outdir}/tmp/{primer_fa_name}
{self.params.blastn} -num_threads {self.params.threads} \
    -query {self.outdir}/{self.prefix}.fasta \
    -db {self.outdir}/tmp/{primer_fa_name} \
    -task blastn-short -word_size {self.params.blastn_wordsize} -evalue {self.params.blastn_evalue} \
    -out {self.outdir}/{self.prefix}.fasta.out \
    -outfmt '{self.params.blastn_outfmt}'
"""
        logging.info(cml)
        run(cml, shell=True)

    def get_dict_primer(self):
        #blastn outfmt
        dic = {fn[1]:fn[0] for fn in enumerate(self.params.blastn_outfmt.split(' ')[1:])}
        Head = namedtuple('Head', dic.keys())
        head = Head(**dic)
        # read-primer字典, {qseqid:[[sseqid,..],['qstart-qend,sstart-send',..]], ...}
        dict_read_primer = defaultdict(list)
        with open(f'{self.outdir}/{self.prefix}.fasta.out') as f:
            for line in f:
                llst = line.strip().split('\t')
                qid, sid = llst[head.qseqid], llst[head.sseqid]
                info = f'{llst[head.qstart]}-{llst[head.qend]},{llst[head.sstart]}-{llst[head.send]}'
                dict_read_primer.setdefault(qid, [[],[]])
                dict_read_primer[qid][0].append(sid)
                dict_read_primer[qid][1].append(info)
        #primer字典, {'primer1-primer2..': [[info1;info2-..],[rid1;rid2-..]}
        dict_primer = defaultdict(list)
        for qid in dict_read_primer:
            if len(dict_read_primer[qid][0]) < 2: continue
            sid = '-'.join(dict_read_primer[qid][0])
            info = ';'.join(dict_read_primer[qid][1])
            dict_primer.setdefault(sid, [[],[]])
            dict_primer[sid][1].append(qid)
            dict_primer[sid][0].append(info)
        return dict_primer
    
    def count_primer(self):
        """统计引物"""
        dict_primer = self.get_dict_primer()
        #输出 二聚体组 位置信息 read名
        with open(f'{self.outdir}/{self.prefix}.primer_groups.tsv', 'wt', encoding='utf-8', newline='') as g:
            g.write('group\tposition\treadid\n')
            for sid in dict_primer:
                for idx in range(len(dict_primer[sid][0])):
                    g.write(f'{sid}\t{dict_primer[sid][0][idx]}\t{dict_primer[sid][1][idx]}\n')
        #二聚体组 位置信息 绘图
        df = pd.read_table(f'{self.outdir}/{self.prefix}.primer_groups.tsv')
        #输出引物组计数文件
        self.df_count = df.iloc[:,:2].groupby('group').count().sort_values(by='position', ascending=False)\
                        .rename(columns={'position':'count'})
        self.df_count.to_csv(f'{self.outdir}/{self.prefix}.primer_groups_count.tsv', sep='\t', encoding='utf-8')
        #引物组位置计数文件
        self.df_plot2 = df.groupby(['group','position']).count().sort_values(by='readid', ascending=False).reset_index()\
                        .rename(columns={'readid':'count'})
        self.df_plot2.to_csv(f'{self.outdir}/{self.prefix}.primer_groups_position_count.tsv', sep='\t', encoding='utf-8', index=False)

    def myplot(self):
        """
        统计图
        1. 整体引物结果的饼图和柱形图
        2. 每个引物组的position柱形图
        """
        #左饼图,右柱形图
        df_count_others = pd.DataFrame([[self.df_count.iloc[9:,]['count'].sum()]], index=['Others'], columns=['count'])
        df_count_plot = pd.concat([df_count_others, self.df_count.iloc[:9,].sort_values(by='count')])
        fig = plt.figure(figsize=(20, 10), facecolor='white')
        fig.suptitle(self.prefix)
        ax1 = fig.add_subplot(121)
        # labels = [df_count_plot.index[0]] + [None]*3 + list(df_count_plot.index[4:])
        ax1.pie(df_count_plot['count'], rotatelabels=True, labels=df_count_plot.index, 
                textprops={'color':'white'}, colors=plt.get_cmap('Set3')(range(10)),
                wedgeprops={'linewidth':1,'edgecolor':'white'})
        ax1.legend(loc='upper left')
        ax2 = fig.add_subplot(122)
        ax2.barh(df_count_plot.index, df_count_plot['count'].values, color=plt.get_cmap('Set3')(range(10)))
        ax2.set_xlabel('count')
        ax2.spines['top'].set_color('none')
        ax2.spines['right'].set_color('none')
        fig.savefig(f'{self.outdir}/{self.prefix}.dimer_analysis.png', dpi=300)
        #positino柱状图
        self.df_plot2.rename(columns={'readid':'count'}, inplace=True)
        dir_bar = Path(self.outdir).joinpath(f'{self.prefix}_barplots')
        dir_bar.mkdir(exist_ok=True)
        for group in df_count_plot[1:].index:
            self.df_plot2[self.df_plot2['group']==group][:20].plot.bar(x='position', y='count', title=group, fontsize=6)
            plt.subplots_adjust(bottom=0.4)
            plt.savefig(f'{dir_bar}/{group}.png',dpi=200)

    def execute(self):
        self.seqtk_blastn()
        self.count_primer()
        self.myplot()
        run(f'rm -r {self.outdir}/tmp', shell=True)


class DimerAssessPE(DimerAssess):
    def __init__(self, fastq1, fastq2, primer_fa, outdir, prefix):
        """[230316] 双端测序去引物"""
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.primer_fa = primer_fa
        self.outdir = outdir
        self.prefix = prefix
        #Others
        self.params = config()
        Path(self.outdir).joinpath('tmp').mkdir(parents=True, exist_ok=True)

    def pandaseq_seqtk_blastn(self):
        """
        1. seqtk fastq转fasta; 
        2. blast引物序列建库, blast比对reads和引物
        """
        primer_fa_name = Path(self.primer_fa).name
        cml = f"""
#合并双端
{self.params.pandaseq} -A flash -T {self.params.threads} -t 0.01 \
    -f {self.fastq1} -r {self.fastq2} \
    -w {self.outdir}/tmp/{self.prefix}.pair.fasta -u {self.outdir}/tmp/{self.prefix}.unpair.fasta \
    -g {self.outdir}/tmp/{self.prefix}.pandaseq.log
cat {self.outdir}/tmp/{self.prefix}.pair.fasta {self.outdir}/tmp/{self.prefix}.unpair.fasta > {self.outdir}/{self.prefix}.fasta
#比对reads到primer
cp {self.primer_fa} {self.outdir}/tmp/{primer_fa_name}
{self.params.makeblastdb} -in {self.outdir}/tmp/{primer_fa_name} -dbtype nucl -out {self.outdir}/tmp/{primer_fa_name}
{self.params.blastn} -num_threads {self.params.threads} \
    -query {self.outdir}/{self.prefix}.fasta \
    -db {self.outdir}/tmp/{primer_fa_name} \
    -task blastn-short -word_size {self.params.blastn_wordsize} -evalue {self.params.blastn_evalue} \
    -out {self.outdir}/{self.prefix}.fasta.out \
    -outfmt '{self.params.blastn_outfmt}'
"""
        logging.info(cml)
        run(cml, shell=True)

    def get_dict_primer(self):
        #blastn outfmt
        dic = {fn[1]:fn[0] for fn in enumerate(self.params.blastn_outfmt.split(' ')[1:])}
        Head = namedtuple('Head', dic.keys())
        head = Head(**dic)
        # read-primer字典, {qseqid:[[sseqid,..],['qstart-qend,sstart-send',..]], ...}
        dict_read_primer = defaultdict(list)
        with open(f'{self.outdir}/{self.prefix}.fasta.out') as f:
            for line in f:
                llst = line.strip().split('\t')
                qid, sid = llst[head.qseqid], llst[head.sseqid]
                info = f'{llst[head.qstart]}-{llst[head.qend]},{llst[head.sstart]}-{llst[head.send]}'
                dict_read_primer.setdefault(qid, [[],[]])
                dict_read_primer[qid][0].append(sid)
                dict_read_primer[qid][1].append(info)
        #primer字典, {'primer1-primer2..': [[info1;info2-..],[rid1;rid2-..]}
        dict_primer = defaultdict(list)
        for qid in dict_read_primer:
            if len(dict_read_primer[qid][0]) == 1: #只有一个引物, 那就是(理论)没有引物二聚体
                continue
            elif (len(dict_read_primer[qid][0]) == 2): 
                prmr1, prmr2 = dict_read_primer[qid][0]
                #[230316]匹配两个引物, 是一对(不能是同一个引物). 一对引物的命名规则: S0033NT2F1-S0033NT2R1
                re.sub('F|R','',prmr1)
                if (re.sub('F|R','',prmr1) == re.sub('F|R','',prmr2)) and (len(list(set([prmr1, prmr2]))) != 1):
                    continue
            sid = '-'.join(dict_read_primer[qid][0])
            info = ';'.join(dict_read_primer[qid][1])
            dict_primer.setdefault(sid, [[],[]])
            dict_primer[sid][1].append(qid)
            dict_primer[sid][0].append(info)
        return dict_primer
    
    def execute(self):
        self.pandaseq_seqtk_blastn()
        self.count_primer()
        self.myplot()
        run(f'rm -r {self.outdir}/tmp', shell=True)
