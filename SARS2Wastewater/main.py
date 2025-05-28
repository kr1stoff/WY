#!/usr/bin/env python

# @CreateTime       : 2023/3/28
# @Author           : mengxf
# @version          : v2.0
# @LastModified     : 2023/4/25
# @description      : 新冠污水监测流程

import click
import logging
# custom
from lib import mypreprocess
from lib import myanalysis
from lib import upload_report

# 运行日志
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('-i', '--inyaml', type=click.Path(exists=True), required=True, help='输入YAML文件, SOLAR生成.')
@click.option('-o', '--analysis', default='.', show_default=True,
              help='结果分析目录, 在该目录下层新建文件夹保存结果. e.g. {analysis}/{task_id}')
@click.option('-p', '--platform', default='ILMN', show_default=True,
              type=click.Choice(['ILMN', 'ONT']), help='测序平台, 支持ILMN,ONT.')
def main(inyaml, analysis, platform):
    """新冠污水监测流程."""
    # 前处理
    prprcs = mypreprocess.PreProcess(inyaml, analysis, platform)
    prprcs.parse_inyaml()
    prprcs.get_params()
    prprcs.make_directories()
    prprcs.link_fastq()
    prprcs.generate_snakemake_config()
    # 分析 ILMN/ONT
    apipe = myanalysis.Analysis(prprcs)
    apipe.run_snakemake_pipe()
    apipe.handle_user_predict_lineages()
    # 上传&报告 ILMN/ONT
    urpipe = upload_report.UploadReport(prprcs)
    if platform == 'ILMN':
        urpipe.upload()
    else:
        urpipe.uploadONT()
    urpipe.report()
    urpipe.zip_upload()


if __name__ == '__main__':
    main()
