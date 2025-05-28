#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/5/23 19:53
# @Last Modified by:   Ming
# @Last Modified time: 2022/5/23 19:53
import logging
import os.path
import sys
from pathlib import Path
from subprocess import run

import click
import pandas as pd
import yaml

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#### Some Setting
PATH = Path(__file__).absolute().parent
LIBPATH = PATH.joinpath("lib")
BINPATH = PATH.joinpath("bin")
CONFIGPATH = PATH.joinpath("config")
from lib.mlst import MLST
from lib.cgmlst import cgMLST
from lib.wgmlst import wgMLST
from lib.ptyping import PTyping
from lib.cgsnp import cgSNP
from lib.wgsnp import wgSNP
from lib.utils import Mutex


#### Some Function
def report(species, d_report, name):
    """
    报告生成与打包
    """
    logger.info("Generate the report and package")
    d_name = Path(d_report)
    cmd = f"""set -e
perl {BINPATH}/report.pl -species {species} -outdir {d_name.absolute()}
cd {d_name.parent}
zip -r {name}.zip {name} >/dev/null 2>&1
"""
    logger.debug(cmd)
    run(cmd, shell=True)


def get_type_info(f_content):
    """
    Get the type information from content file
    """
    df = pd.read_csv(f_content, sep="\t")
    return df


def get_xsnp_ref(pathogen):
    """
    Get the xSNP reference name

    :param pathogen: The pathogen name
    """
    f_database = CONFIGPATH.joinpath("database.yml")
    info_database = yaml.safe_load(open(f_database, 'r').read())
    return Path(info_database["typetrace"]).joinpath(f"xSNP/{pathogen}/ref.fna")


#### Main
f_content = CONFIGPATH.joinpath("content.tsv")
info_content = get_type_info(f_content)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-c", "--config",
              type=click.Path(),
              cls=Mutex,
              not_required_if=["samples", "fastas", "db"],
              help="The VENUS yaml config file")
@click.option("--samples",
              type=click.Path(),
              cls=Mutex,
              not_required_if=["config"],
              help="The sample name sep by ,")
@click.option("--fastas",
              type=click.Path(),
              cls=Mutex,
              not_required_if=["config"],
              help="The fasta genome sequence sep by ,")
@click.option('--pathogen',
              cls=Mutex,
              not_required_if=["config"],
              type=click.Choice(list(info_content["pathogen_name"])),
              help="The pathogen name to use")
@click.option('-o', '--out',
              cls=Mutex,
              not_required_if=["config"],
              type=click.Path(),
              help="The out put dir")
@click.option('--upload',
              type=click.Path(),
              default="",
              show_default=True,
              help="The AnalysisResults + taskid dir")
@click.option('--mlst',
              cls=Mutex,
              not_required_if=["config"],
              type=bool,
              default=True,
              help="Wheather run the mlst analysis")
@click.option('--cgmlst',
              cls=Mutex,
              not_required_if=["config"],
              type=bool,
              default=True,
              help="Wheather run the cgmlst analysis")
@click.option('--wgmlst',
              cls=Mutex,
              not_required_if=["config"],
              type=bool,
              default=True,
              help="Wheather run the wgmlst analysis")
@click.option('--cgsnp',
              cls=Mutex,
              not_required_if=["config"],
              type=str,
              default=None,
              help="Wheather run the cgsnp analysis and the database to use")
@click.option('--wgsnp',
              cls=Mutex,
              not_required_if=["config"],
              type=str,
              default=None,
              help="Wheather run the wgsnp analysis and the database to use")
@click.option('--serotype',
              cls=Mutex,
              not_required_if=["config"],
              type=bool,
              default=True,
              help="Wheather run the serotype analysis")
@click.option('--other',
              cls=Mutex,
              not_required_if=["config"],
              type=str,
              default=None,
              help="Wheather run the other typing analysis and the database to use, sep by ;")
@click.option("--debug",
              is_flag=True,
              show_default=True,
              help="Whether turn on the debug mode")
def main(config, samples, fastas, pathogen, out, upload, mlst, cgmlst, wgmlst, cgsnp, wgsnp, serotype, other, debug):
    """
    一体机分型溯源流程

    如提供一体机生成的配置文件，则不需要其它参数\n

    TODO: cgSNP,wgSNP参考的选择
    """
    if debug:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # 获取数据库路径
    f_database = CONFIGPATH.joinpath("database.yml")
    info_database = yaml.safe_load(open(f_database, 'r').read())

    # 配置信息获取
    if config:
        # 来自一体机
        f_venus = Path(config).absolute()
        info_venus = yaml.safe_load(open(f_venus, 'r').read())
        list_samples = list(info_venus["samples"].keys())
        list_fastas = [info_venus["samples"][i]["fasta"] for i in list_samples]
        out = Path(f"/sdbb/Earth/Analysis/{info_venus['task_id']}").absolute()
        upload = Path(f"/sdbb/Earth/AnalysisResults/{info_venus['task_id']}").absolute() if upload == "" else None
        pathogen_id = int(info_venus["genome_typing_database"])
        pathogen_info = info_content[info_content["pathogen_id"] == pathogen_id].reset_index()
        logger.debug(pathogen_info)
        pathogen_name = pathogen_info.loc[0, "pathogen_name"]
        scientific_name = pathogen_info.loc[0, "scientific_name"]
        mlst = True if info_venus["mlst"] else False
        cgmlst = True if info_venus["cgmlst"] else False
        wgmlst = True if info_venus["wgmlst"] else False
        cgsnp = True if info_venus["cgsnp"] else False
        wgsnp = True if info_venus["wgsnp"] else False
        serotype = True if info_venus["serotype"] else False
        other = info_venus["other_type"] if info_venus["other_type"] else False
    else:
        pathogen_name = pathogen
        pathogen_info = info_content[info_content["pathogen_name"] == pathogen]
        list_samples = samples.strip().split(',')
        list_fastas = fastas.strip().split(',')
        out = Path(out).absolute()

    # 项目分型结果统计
    types = ["mlst", 'cgmlst', 'wgmlst', 'serotype', 'other_type']
    dict_type_result = {sample: {i: "" for i in types} for sample in list_samples}

    # 目录生成
    d_upload = out.joinpath("Upload")
    out.mkdir(exist_ok=True, parents=True)
    d_upload.mkdir(exist_ok=True, parents=True)

    # Start to analysis
    ## xSNP分析，需要多个样本
    ### cgSNP
    if cgsnp and pathogen_info.loc[0, "cgsnp"] and len(list_samples) >= 3:
        logger.info("cgSNP analysis")
        d_out = out.joinpath("cgSNP")
        d_out.mkdir(exist_ok=True, parents=True)
        f_sample = d_out.joinpath("SampleSheet.tsv")
        with open(f_sample, 'w') as OUT:
            for i, j in zip(list_samples, list_fastas):
                print(*[i, j], sep="\t", file=OUT)
        f_ref = get_xsnp_ref(scientific_name.replace(" ", "_"))
        f_mask = f_ref.parent.joinpath("mask.bed")
        job_cgsnp = cgSNP(f_sample, f_ref, f_mask, d_out)
        job_cgsnp.run()
        for name in list_samples:
            d_report = d_upload.joinpath(name)
            d_report.mkdir(exist_ok=True, parents=True)
            job_cgsnp.to_upload(d_report)
    ### wgSNP
    if wgsnp and pathogen_info.loc[0, "wgsnp"] and len(list_samples) >= 3:
        logger.info("wgSNP analysis")
        d_out = out.joinpath("wgSNP")
        d_out.mkdir(exist_ok=True, parents=True)
        f_sample = d_out.joinpath("SampleSheet.tsv")
        with open(f_sample, 'w') as OUT:
            for i, j in zip(list_samples, list_fastas):
                print(*[i, j], sep="\t", file=OUT)
        f_ref = get_xsnp_ref(scientific_name.replace(" ", "_"))
        job_wgsnp = wgSNP(f_sample, f_ref, d_out)
        job_wgsnp.run()
        for name in list_samples:
            d_report = d_upload.joinpath(name)
            d_report.mkdir(exist_ok=True, parents=True)
            job_wgsnp.to_upload(d_report)

    ## 分样本分析的内容
    for name, fasta in zip(list_samples, list_fastas):
        if not os.path.exists(fasta):
            logger.error(sys.exit(f"{fasta} doesn't exists"))
        d_out = out.joinpath(name)
        d_out.mkdir(exist_ok=True, parents=True)
        d_report = d_upload.joinpath(name)
        d_report.mkdir(exist_ok=True, parents=True)
        # MLST
        if mlst and pathogen_info.loc[0, "mlst"]:
            logger.info("MLST analysis")
            scheme = pathogen_info.loc[0, "mlst"]
            job_mlst = MLST(fasta, scheme, d_out.joinpath("mlst"))
            job_mlst.run()
            job_mlst.to_upload(d_report)
            dict_type_result[name]["mlst"] = job_mlst.result
        # cgMLST
        if cgmlst and pathogen_info.loc[0, "cgmlst"]:
            logger.info("cgMLST analysis")
            scheme = pathogen_info.loc[0, "cgmlst"]
            job_cgmlst = cgMLST(fasta, scheme, d_out.joinpath("cgmlst"))
            job_cgmlst.run()
            job_cgmlst.to_upload(d_report)
            dict_type_result[name]["cgmlst"] = job_cgmlst.result
        # wgMLST
        if wgmlst and pathogen_info.loc[0, "wgmlst"]:
            logger.info("wgMLST analysis")
            scheme = pathogen_info.loc[0, "wgmlst"]
            job_wgmlst = wgMLST(fasta, scheme, d_out.joinpath("wgmlst"))
            job_wgmlst.run()
            job_wgmlst.to_upload(d_report)
            dict_type_result[name]["wgmlst"] = job_wgmlst.result
        # SeroType
        if serotype and pathogen_info.loc[0, "serotype"]:
            logger.info("SeroType analysis")
            job_typing = PTyping(fasta, pathogen_name, d_out.joinpath("serotype"))
            job_typing.run()
            job_typing.to_upload(d_report)
            dict_type_result[name]["serotype"] = job_typing.result
        # Other
        if other and pathogen_info.loc[0, "other_type"]:
            dict_type_result[name]["other_type"] = []
            for i in other.strip().split(';'):
                if i in pathogen_info.loc[0, "other_type"]:
                    if i == "MLVA":
                        logger.info("MLVA analysis")
                        from lib.mlva import MLVA
                        scientific_name = pathogen_info.loc[0, "scientific_name"]
                        pathogen_name.replace(" ", "_")
                        f_primer = f'{info_database["mlva"]}/{scientific_name.replace(" ", "_")}/primers.txt'
                        job_mlva = MLVA(Path(fasta), Path(f_primer), d_out.joinpath("mlva"))
                        job_mlva.run()
                        job_mlva.to_upload(d_report)
                    elif i == "SpaType":
                        logger.info("SpaType analysis")
                        d_spa = d_out.joinpath("spa")
                        d_spa.mkdir(exist_ok=True, parents=True)
                        d_spa_ref = f'{info_database["other"]}/Staphylococcus_aureus/Spa'
                        from lib.spa import Spa
                        job_spa = Spa(fasta, f"{d_spa_ref}/primers.fna", f"{d_spa_ref}/spa_db.csv",
                                      f"{d_spa_ref}/spatypes.txt", d_spa)
                        job_spa.run()
                        job_spa.to_upload(d_report)
                        dict_type_result[name]["other_type"].append(f"SpaType:{job_spa.result}")
                    elif i == "canSNP":
                        logger.info("canSNP analysis")
                        d_cansnp = d_out.joinpath("cansnp")
                        d_cansnp.mkdir(exist_ok=True, parents=True)
                        d_cansnp_ref = f"{info_database['other']}/Bacillus_anthracis/CanSNP"
                        from lib.cansnp import CanSNP
                        job_cansnp = CanSNP(fasta, d_cansnp_ref, f"{d_cansnp_ref}/ref.db", d_cansnp)
                        job_cansnp.run()
                        job_cansnp.to_upload(d_report)
                        dict_type_result[name]["other_type"].append(f"canSNP:{job_cansnp.result}")
        # 报告生成
        report(pathogen_name, d_report, name)
    if upload:
        logger.info(f"Copy the result to {upload}")
        cmd = f"cp -rf {d_upload}/* {upload}"
        run(cmd, shell=True)

    # 分型结果统计
    f_type_stat = f"{d_upload}/type_stat.txt"
    logger.info(f"Generate the type info of this project to {f_type_stat}")
    with open(f_type_stat, 'w') as OUT:
        header = ["Sample", *types, "all_class_type"]
        print(*header, sep="\t", file=OUT)
        for sample in list_samples:
            dict_type_result[sample]["other_type"] = ';'.join(dict_type_result[sample]["other_type"])
            info = [sample]
            type_method = []
            for i in types:
                info.append(dict_type_result[sample][i])
                if i != "other_type" and dict_type_result[sample][i] != "":
                    type_method.append(f"{i}:{dict_type_result[sample][i]}")
            if dict_type_result[sample]["other_type"] != "":
                type_method.append(dict_type_result[sample]["other_type"])
            type_method = ';'.join(type_method)
            info.append(type_method)
            print(*info, sep="\t", file=OUT)


if __name__ == "__main__":
    main()
