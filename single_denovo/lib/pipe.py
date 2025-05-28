#!/home/earthtest/miniconda3/bin/python
# -*- coding: utf-8 -*-
# @Author: zhuzi
# @Date:   2023/5/10 
import os
import sys
import logging
import yaml

logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s: %(message)s',datefmt='%Y-%m-%d %H:%M',stream=sys.stdout)
sys.path.append(os.path.dirname(__file__))

class Pipe():
    def __init__(self,solar_yaml,id,ref):
        # === 读取配置文件 =========================================================
        from lib.common import config
        self.params = config()
        self.Script = os.path.dirname(os.path.dirname(__file__))   # single denovo 主目录
        self.Bin = os.path.join(self.Script,'bin')

        # === 参数 =========================================================
        self.id = id
        self.solar_yaml = os.path.abspath(solar_yaml)
        self.ref = os.path.abspath(ref) if ref != 'None' else False


    def read_solar_yaml(self):
        """读取solar yaml"""
        with open(self.solar_yaml,'r',encoding='utf-8') as file:
            self.task_yaml = yaml.load(file.read(), Loader=yaml.FullLoader)


    def get_path(self):
        """获取分析目录/上传目录 ,创建初始目录"""
        task_id = str(self.task_yaml['task_id'])
        self.out = os.path.join(self.params.Analysis, task_id, self.id)
        self.log = os.path.join(self.params.Analysis, task_id, self.id, 'log')
        self.p_sh = os.path.join(self.params.Analysis, task_id, self.id, '0.shell')
        self.upload = os.path.join(self.params.Upload, task_id)
        os.makedirs(self.p_sh, exist_ok=True)
        os.makedirs(self.upload, exist_ok=True)
        os.makedirs(self.log, exist_ok=True)


    def get_mode(self):
        """获取测序模式 PE/SE"""
        self.PE_flag = True if self.task_yaml['pe_or_se'] == '双端' else False


    def get_fq(self):
        """获取fq路径"""
        if self.PE_flag:
            self.fq1 = self.task_yaml['samples'][self.id]['sequencing_data1']
            self.fq2 = self.task_yaml['samples'][self.id]['sequencing_data2']
        else:
            self.fq1 ,self.fq2 = self.task_yaml['samples'][self.id]['sequencing_data1'] ,''


    def get_pathogen_type(self):
        """获取样本类型： 细菌/真菌/病毒"""
        if self.task_yaml['pathogen_type'] == '细菌':
            self.kind = 'bacteria'
        elif self.task_yaml['pathogen_type'] == '真菌':
            self.kind = 'fungi'
        elif self.task_yaml['pathogen_type'] == '病毒':
            self.kind = 'virus'


    def get_limit_opt(self):
        """获取暂定参数，去污染 / 病毒无参组装，同一个参数"""
        if self.task_yaml['sars_cov_2_database'] == "是":
            self.if_pollute = True
            self.if_denovo = True
        else:
            self.if_pollute = False
            self.if_denovo = False


    def get_var(self):
        """整合关键性变量"""
        # ============================================================================
        self.p_qc = os.path.join(self.out,'1.qc')
        self.p_fastp = os.path.join(self.out,'1.qc','fastp')
        self.p_fastqc = os.path.join(self.out,'1.qc','fastqc')
        self.p_filter_fastqc = os.path.join(self.out,'1.qc','filter_fastqc')

        # ============================================================================
        self.p_ass = os.path.join(self.out,'2.ass')
        self.p_deal = os.path.join(self.out,'2.ass','deal_reads')
        self.p_spades = os.path.join(self.out,'2.ass','spades')
        self.p_quast = os.path.join(self.out,'2.ass','quast')
        self.p_depth = os.path.join(self.out,'2.ass','depth')
        self.p_checkm = os.path.join(self.out,'2.ass','checkm')
        self.p_ref = os.path.join(self.out,'2.ass','ref')

        # ============================================================================
        self.p_pre = os.path.join(self.out,'3.pre')
        self.p_predict = os.path.join(self.out,'3.pre','predict')
        self.p_phispy = os.path.join(self.out,'3.pre','phispy')
        self.p_island = os.path.join(self.out,'3.pre','island')

        # ============================================================================
        self.p_anno = os.path.join(self.out,'4.anno')
        self.p_drug = os.path.join(self.out,'4.anno','drug')
        self.p_VFDB = os.path.join(self.out,'4.anno','VFDB')
        # self.p_pfam = os.path.join(self.out,'4.anno','pfam')
        self.p_swissprot = os.path.join(self.out,'4.anno','swissprot')
        self.p_eggmapper = os.path.join(self.out,'4.anno','eggmapper')
        self.p_CAZy = os.path.join(self.out,'4.anno','CAZy')

        # ============================================================================
        self.p_report = os.path.join(self.out,self.id)

        # == 各模块需要的交叉变量 =====================================================
        if self.PE_flag:
            # clean reads
            self.clean_fq1 = os.path.join(self.p_fastp,f"{self.id}_1.clean.fq.gz")
            self.clean_fq2 = os.path.join(self.p_fastp,f"{self.id}_2.clean.fq.gz")
            self.clean_fq = self.clean_fq1 + ',' +self.clean_fq2
            # 去污染序列
            self.depollute_fq1 = os.path.join(self.p_deal,'depollute_1.fq')
            self.depollute_fq2 = os.path.join(self.p_deal,'depollute_2.fq')
            self.depollute_fq = self.depollute_fq1 + ',' +self.depollute_fq2
            # 截取序列
            self.cutfq1 = os.path.join(self.p_deal,'cut_1.fq')
            self.cutfq2 = os.path.join(self.p_deal,'cut_2.fq')
        
        else:
            self.clean_fq1 = os.path.join(self.p_fastp,f"{self.id}_1.clean.fq.gz")
            self.clean_fq2 = ''
            self.clean_fq = self.clean_fq1

            self.depollute_fq1 = os.path.join(self.p_deal,'depollute_1.fq')
            self.depollute_fq2 = ''
            self.depollute_fq = self.depollute_fq1

            self.cutfq1 = os.path.join(self.p_deal,'cut_1.fq')
            self.cutfq2 = ''
        
        # 过滤前/过滤后，组装序列
        self.raw_contig = os.path.join(self.p_spades,'scaffolds.fasta')
        self.filter_contig = os.path.join(self.p_spades,f'{self.id}_scaffolds.fasta')
        # 组装序列比对原始fq的bam
        self.sort_bam = os.path.join(self.p_depth,'sort.bam')
        # 预测核酸fa，蛋白faa ，用于后续注释使用
        self.predict_fa = os.path.join(self.p_predict,'predict.ffn')
        self.predict_faa = os.path.join(self.p_predict,'predict.faa')
        # 注释需要的深度文件，统计了各个序列的深度/覆盖度
        self.merge_file = os.path.join(self.p_predict,'bed_bam.merge_depth_cov.txt')



    def dic2yaml(self):
        """方便单模块运行"""
        pipeline_yaml = os.path.join(self.out,'pipeline.yaml')
        # 变量字典
        var_dic = {
            'ref': self.ref,
            'PE_flag': self.PE_flag,
            'if_pollute': self.if_pollute,
            'if_denovo': self.if_denovo,
            'id': self.id,
            'kind': self.kind,
            'out': self.out,
            'upload': self.upload,
            'p_sh': self.p_sh,
            'log': self.log,

            'p_qc': self.p_qc,
            'p_fastp': self.p_fastp,
            'p_fastqc': self.p_fastqc,
            'p_filter_fastqc': self.p_filter_fastqc,

            'p_ass': self.p_ass,
            'p_deal': self.p_deal,
            'p_spades': self.p_spades,
            'p_quast': self.p_quast,
            'p_checkm': self.p_checkm,
            'p_depth': self.p_depth,
            'p_ref': self.p_ref,

            'p_pre': self.p_pre,
            'p_predict': self.p_predict,
            'p_island': self.p_island,
            'p_phispy': self.p_phispy,

            'p_anno': self.p_anno,
            'p_drug': self.p_drug,
            'p_VFDB': self.p_VFDB,
            'p_eggmapper': self.p_eggmapper,
            'p_swissprot': self.p_swissprot,
            'p_CAZy': self.p_CAZy,

            'p_report': self.p_report,

            'fq1': self.fq1,
            'fq2': self.fq2,
            'clean_fq1': self.clean_fq1,
            'clean_fq2': self.clean_fq2,
            'depollute_fq1': self.depollute_fq1,
            'depollute_fq2': self.depollute_fq2,
            'cutfq1': self.cutfq1,
            'cutfq2': self.cutfq2,
            'raw_contig': self.raw_contig,
            'filter_contig': self.filter_contig,
            'sort_bam': self.sort_bam,
            'predict_fa': self.predict_fa,
            'predict_faa': self.predict_faa,
            'merge_file': self.merge_file

        }
        # 将上述字典变量，写入样本的分析目录
        if not os.path.exists(pipeline_yaml):
            logging.info(f"[Output]: dic===>yaml {pipeline_yaml}")
            with open(pipeline_yaml, 'w' ,encoding='utf-8') as file:
                file.write(yaml.dump(var_dic, allow_unicode=True))


    def run(self):
        """运行本类的函数"""
        self.read_solar_yaml()
        self.get_path()
        self.get_mode()
        self.get_pathogen_type()
        self.get_fq()
        self.get_limit_opt()
        self.get_var()
        self.dic2yaml()

