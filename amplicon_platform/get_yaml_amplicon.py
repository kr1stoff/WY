# -*- coding: utf-8 -*-
import yaml
import click


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-r', '--run', default='nested', type=str, help="运行步骤")
@click.option('-o', '--out', default='./input.yaml', type=click.Path(), help="输出的yaml,默认为./test.yaml")
def main(run, out):

    steps = list(run)
    first_step = steps[0]
    last_step = steps[-1]

    data_with_comments = {}
    data_with_comments['run'] = run

    if (int(first_step) == 1):
        data_with_comments['taxid'] = ''

        data_with_comments['ptype'] = ''
        if (int(last_step) >= 2):
            data_with_comments['resion_type'] = ''
            data_with_comments['gff'] = ''
            # data_with_comments['taxid'] = ''
            data_with_comments['windows'] = '400'
            data_with_comments['step'] = '200'
            data_with_comments['inc_num'] = '20000'
            data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'

        if (int(last_step) >= 3):
            data_with_comments['reference_dir'] = ''
            data_with_comments['name2taxid'] = '/sdbb/bioinfor/lanlei/script/amplicon_platform/src/name2taxid.tsv'
            # data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
            data_with_comments['evalue'] = '800000'
            data_with_comments['base_num'] = '20'
            data_with_comments['chinese'] = ''

        # if (int(last_step) >= 4):
        #     data_with_comments['reference_dir'] = ''

        if (int(last_step) >= 5):
            # data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
            data_with_comments['database_path'] = '/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge'

    if (int(first_step) == 2):
        data_with_comments['reference_dir'] = ''
        data_with_comments['resion_type'] = ''
        data_with_comments['gff'] = ''
        data_with_comments['taxid'] = ''
        data_with_comments['windows'] = '400'
        data_with_comments['step'] = '200'
        data_with_comments['inc_num'] = '20000'
        data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
        if (int(last_step) >= 3):
            # data_with_comments['reference_dir'] = ''
            data_with_comments['name2taxid'] = '/sdbb/bioinfor/lanlei/script/amplicon_platform/src/name2taxid.tsv'
            # data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
            data_with_comments['evalue'] = '800000'
            data_with_comments['base_num'] = '20'
            data_with_comments['chinese'] = ''

        # if (int(last_step) >= 4):
        #     data_with_comments['reference_dir'] = ''

        if (int(last_step) >= 5):
            # data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
            data_with_comments['database_path'] = '/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge'

    if (int(first_step) == 3):
        data_with_comments['bed'] = ''
        data_with_comments['reference_dir'] = ''
        data_with_comments['name2taxid'] = '/sdbb/bioinfor/lanlei/script/amplicon_platform/src/name2taxid.tsv'
        data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
        data_with_comments['evalue'] = '800000'
        data_with_comments['base_num'] = '20'
        data_with_comments['chinese'] = ''

        # if (int(last_step) >= 4):
        #     data_with_comments['reference_dir'] = ''

        if (int(last_step) >= 5):
            # data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
            data_with_comments['database_path'] = '/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge'

    if (int(first_step) == 4):
        data_with_comments['merge_primer'] = ''
        data_with_comments['reference_dir'] = ''

        if (int(last_step) >= 5):
            data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
            data_with_comments['database_path'] = '/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge'

    if (int(first_step) == 5):
        data_with_comments['primer2taxid_table'] = ''
        data_with_comments['nt_database'] = '/sdbb/bioinfor/yehui/nt/nt.fa'
        data_with_comments['database_path'] = '/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge'

    comments = {
        'run': 'Run step, eg:23456\n# 01.download,02.conserved_sequences,03.primer_evaluate,04.degenerate_base,05.primer_reevaluate',
        'taxid': '物种的taxid',
        'ptype': '[Bacteria|Viruses|Fungi]',
        'resion_type': '保守区域范围, [gene|genome]',
        'gff': 'gff文件, 保守区域预测选择基因才需要',
        'windows': '保守区域鉴定时, 滑动窗口大小, 默认400',
        'step': '保守区域鉴定时, 滑动窗口步长, 默认200',
        'inc_num': '保守区域保留的区域数, 默认为20000',
        'nt_database': 'nt库路径, 默认为/sdbb/bioinfor/yehui/nt/nt.fa',
        'reference_dir': '物种基因组数据库路径，一般在/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge下面',
        'bed': '保守区域bed文件',
        'name2taxid': '拉丁名和taxid的对应表, 默认为/sdbb/bioinfor/lanlei/script/amplicon_platform/src/name2taxid.tsv',
        'evalue': '引物比对evalue值, 默认为800000',
        'base_num': 'EXCEL扩增子片段延升的碱基数, 默认为20',
        'chinese': '物种中文名, 输出结果需要',
        'merge_primer': '引物评估后的xls表',
        'primer2taxid_table': '引物和taxid的对应表',
        'database_path': '数据库路径, 默认为/mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge'
    }

    key_order = ['run', 'taxid', 'ptype', 'resion_type', 'gff', 'windows', 'step', 'inc_num', 'nt_database', 'reference_dir',
                 'bed', 'name2taxid', 'evalue', 'base_num', 'chinese', 'merge_primer', 'primer2taxid_table', 'database_path']

    # 构建 YAML 文件内容
    yaml_content = ""
    for key in key_order:
        if key in data_with_comments:
            if key in comments:
                yaml_content += f"# {comments[key]}\n"
            yaml_content += f"{key}: {data_with_comments[key]}\n"

    # 将内容写入 YAML 文件
    with open(out, 'w') as yaml_file:
        yaml_file.write(yaml_content)


if __name__ == '__main__':
    main()
