import pandas as pd
from pathlib import Path
import logging


def get_ref_acc(download_genome_dir, tid):
    """获取参考基因组的accession"""
    ref_genome_name = Path(f'{download_genome_dir}/{tid}/ref_genome_name')
    if (not ref_genome_name.is_file()) or (ref_genome_name.stat().st_size == 0):
        logging.error(f'文件不存在或文件为空: {ref_genome_name} !')
    with open(ref_genome_name) as f:
        acc = '_'.join(f.read().strip().split('_')[:2])
    return acc


def get_remove_assembly_level_accs(download_genome_dir, tid):
    """
    去除assembly_level为Contig/Scaffold后的assembly_accession号,需要保留ref_genome_name
    assembly_info   ::  当前assembly_info.tsv文件
    """
    assembly_info = Path(download_genome_dir).joinpath(f'{tid}/assembly_info.tsv')
    df = pd.read_table(assembly_info, sep='\t', low_memory=False)
    df = df.loc[~df['assembly_level'].isin(['Contig', 'Scaffold'])]
    if df.shape[0] > 2000: # 只保留Complete Genome
        df = df.loc[df['assembly_level'].isin(['Complete Genome'])]
    if df.shape[0] > 2000: # 只保留时间2020后
        df['seq_rel_date'] = pd.to_datetime(df['seq_rel_date'])
        df = df.loc[df['seq_rel_date'].gt(pd.Timestamp('2020-1-1'))]
    accs = list(df['X..assembly_accession'].values)
    #如果参考accession不在筛选后列表里,加进去
    acc = get_ref_acc(download_genome_dir, tid)
    if acc not in accs: accs.append(acc)
    return accs


def acc2txt(accs, txt):
    """accession列表写入文档中"""
    with open(txt, 'wt', newline='', encoding='utf-8') as g:
        for acc in accs:
            g.write(acc + '\n')


def repfunc(work_dir, tid, fna, accessions):
    """软连接并添加accession到列表"""
    if not Path(work_dir).joinpath(f'{tid}/fnas/{fna.name}').exists():
        Path(work_dir).joinpath(f'{tid}/fnas/{fna.name}').symlink_to(fna)
    acc = '_'.join(Path(fna).name.split('_')[:2])
    accessions.append(acc)


def taxid_accession(download_genome_dir, work_dir, tid):
    """获取taxid下的accession号列表"""
    accessions = []
    fnas = list(Path(download_genome_dir).glob(f"{tid}/GC*fna"))
    if len(fnas) <= 2000:
        acc = get_ref_acc(download_genome_dir, tid) #顺便检查一下ref_genome_name
        for fna in fnas:
            repfunc(work_dir, tid, fna, accessions)
    else:
        accs = get_remove_assembly_level_accs(download_genome_dir, tid)
        for acc in accs:
            res = list(Path(download_genome_dir).glob(f'{tid}/{acc}*.fna'))
            if (not res): continue
            fna = res[0]
            repfunc(work_dir, tid, fna, accessions)
    return accessions


def check_completed(work_dir, tid):
    """
    检查,返回布尔值
    1. 已经有'accessions.txt'且不是空文件;
    2. 已经预测过的物种"""
    accessions_file = Path(f'{work_dir}/{tid}/accessions.txt')
    if (accessions_file.is_file()) and (accessions_file.stat().st_size != 0):
        logging.info(f'已经生成accessions_file. {work_dir}/{tid}/accessions.txt')
        return True
    elif Path(f'{work_dir}/{tid}/core_cluster.bed').is_file():
        logging.info(f'已经预测完成. {work_dir}/{tid}/core_cluster.bed')
        return True
    else:
        return False


def accession_list(excel, download_genome_dir, work_dir, skip_completed):
    """
    将过滤后的ACCESSION号写入记录文件中, 软连接各物种fna到work_dir.
    excel               ::  钉钉项目进度表, 表头必有 分类/过滤完成/taxid  
    download_genome_dir ::  梓伦/超博下载基因组的文件目录  
    work_dir            ::  工作目录, prokka/roary/conserve都在这个目录做  
    """
    df = pd.read_excel(excel, sheet_name='Sheet1', dtype=str) #第二批
    taxids = df.loc[(df['病原体类型'].isin(['细菌', '支/衣/立'])) & (df['过滤完成']=='√'), 'TaxonID'].values
    #进度表
    for tid in taxids:
        if skip_completed and check_completed(work_dir, tid): continue
        Path(work_dir).joinpath(f"{tid}/fnas").mkdir(parents=True, exist_ok=True)
        accessions = taxid_accession(download_genome_dir, work_dir, tid)
        acc2txt(accessions, f'{work_dir}/{tid}/accessions.txt')
