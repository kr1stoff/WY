from pathlib import Path
import pandas as pd


def df_schedule_prokka(tid, fna_number, gff_number, download_genome_dir):
    """指定taxid下prokka完成状态"""
    if not Path(download_genome_dir).joinpath(f'{tid}/assembly_info.tsv').exists():
        return 'N'
    if fna_number != gff_number:
        prokka_status = 'N'
    else:
        prokka_status = 'Y'
    return prokka_status


def df_schedule_addition1(df_schedule, work_dir):
    """df_schedule新增列, roary&conserved完成状态"""
    df_schedule['Roary'] = 'N'
    df_schedule['Conserved'] = 'N'
    df_schedule['Bed'] = ''
    for row in df_schedule.iterrows():
        tid = row[1]['taxid']
        if Path(work_dir).joinpath(f'{tid}/roary_out/gene_presence_absence.Rtab').exists():
            df_schedule.loc[row[0], 'Roary'] = 'Y'
        if Path(work_dir).joinpath(f'{tid}/core_cluster.bed').exists():
            df_schedule.loc[row[0], 'Conserved'] = 'Y'
            df_schedule.loc[row[0], 'Bed'] = Path(f'{work_dir}/{tid}/core_cluster.bed').resolve()
    return df_schedule


def predict_schedule(excel, download_genome_dir, work_dir):
    """
    获取预测进度表.
    excel = Path(work_dir).joinpath('materials/221214扩增子数据库下载进度.xlsx')
    """
    df = pd.read_excel(excel, sheet_name='Sheet1', dtype=str) #第二批
    taxids = df.loc[(df['病原体类型'].isin(['细菌', '支/衣/立'])) & (df['过滤完成']=='√'), 'TaxonID'].values
    #进度表
    #prokka
    records = []
    for tid in taxids:
        raw_fna_number = len(list(Path(download_genome_dir).glob(f'{tid}/GC*.fna')))
        gff_number = len(list(Path(work_dir).glob(f'{tid}/GC*.gff')))
        fna_number = len(list(Path(work_dir).glob(f'{tid}/fnas/GC*.fna')))
        d_status = 'Y' if Path(download_genome_dir).joinpath(f'{tid}/assembly_info.tsv').exists() else 'N'
        p_status = df_schedule_prokka(tid, fna_number, gff_number, download_genome_dir)
        records.append([tid, raw_fna_number, fna_number, gff_number, d_status, p_status])
    df_schedule = pd.DataFrame(records, columns=['taxid', 'raw_fna_number', 'fna_number', 'gff_number', 'Download', 'Prokka'])
    df_schedule = df_schedule_addition1(df_schedule, work_dir)
    df_schedule.to_csv(f'{work_dir}/predict_schedule.tsv', sep='\t', encoding='utf-8', index=False)
