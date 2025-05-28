import os
import time
from subprocess import run

import pandas as pd
import pymysql
import yaml


#############################
# mysql 操作
# db.cursor()
# execute() 执行
# fetchone() 方法获取单条数据
# fetchall() 方法获取多条数据
#############################

def return2mysql(list1,
                 cursor,
                 task_id,
                 zip_pattern=r"/sdbb/Earth/AnalysisResults/%s/%s.zip",
                 html_pattern=r"http://172.16.0.18:8080/file/weiyuan/AnalysisResults/%s/%s/index.html",
                 fa_pattern=""):
    """
    输入: 1.样本名列表 [A1,A2,A3,...]
          2.MySQL游标对象
          3.任务ID
          3.单样本结果压缩包路径通配,第一个%s是任务ID,第二个%s是样本ID
          4.单样本网页结果路径通配
          5.单样本结果FASTA序列路径,默认""无FASTA结果
    结果更新到MySQL数据库上,对应任务ID和样本ID的返回位置
    """
    for sample_id in list1:
        # zip
        return_file_path = zip_pattern % (task_id, sample_id)
        sql_path = "update tb_task_sample set return_file_path='%s' where task_id='%s' and sample_number='%s'" % (
            return_file_path, task_id, sample_id)
        print(sql_path)
        cursor.execute(sql_path)
        # index html
        return_html_url = html_pattern % (task_id, sample_id)
        sql_url = "update tb_task_sample set return_html_url='%s' where task_id='%s' and sample_number='%s'" % (
            return_html_url, task_id, sample_id)
        print(sql_url)
        cursor.execute(sql_url)
        # sequencing fasta
        if fa_pattern != "":
            # /sdbb/Earth/Analysis/2333/XG0120/4.consensus/XG0120.consensus.fa
            # /sdbb/Earth/Analysis/%s/%s/4.consensus/consensus.fa
            return_fasta = fa_pattern % (task_id, sample_id)
            sql_fasta = "update tb_task_sample set return_fasta='%s' where task_id='%s' and sample_number='%s'" % (
                return_fasta, task_id, sample_id)
            cursor.execute(sql_fasta)


# Add by renchaobo
def linkdata(info_sample, d_rawdata):
    """
    Link the rawdata to the analysis dir and rename it to the standard format

    :param info_sample: The sample info table
    :param d_rawdata: The rawdata dir in the pipeline
    """
    samples = []
    for record in info_sample:
        sample_name = record[7]
        sample_type = record[21]
        if sample_type == "PE150":
            cmd = f"""ln -sf {record[23]} {d_rawdata}/{sample_name}_1.fq.gz
ln -sf {record[24]} {d_rawdata}/{sample_name}_2.fq.gz"""
            run(cmd, shell=True)
        elif sample_type == "SE50":
            cmd = f"ln -sf {record[23]} {d_rawdata}/{sample_name}_1.fq.gz"
            run(cmd, shell=True)
        else:
            raise ValueError(f"UNSUPPORT SAMPLETYPE: {sample_type}")
        samples.append(sample_name)
    return samples


def get_group_info(info_sample):
    """
    Get the group info

    :param info_sample: The sample info table
    """
    res = {}
    for record in info_sample:
        sample_name = record[7]
        group_name = record[19]
        if group_name:
            real_group_name = group_name.strip().split('(')[1].strip(')')
            res.setdefault(real_group_name, [])
            res[real_group_name].append(sample_name)
    if len(res) == 0:
        return None
    else:
        return res


def get_diff_info(info_sample):
    """
    Get the group diff info

    :param info_sample: The sample info table
    """
    res = {}
    for record in info_sample:
        info_compare = record[38]
        if info_compare:
            res[tuple(info_compare)] = 1

    if len(res) > 0:
        return [list(i) for i in res.keys()]
    else:
        return None


def config4amplicon(info_task, info_sample, d_analysis):
    """
    Generate YAML config file for 16S/18S/ITS

    :param info_task: The task info get from table tb_analyse_task
    :param info_sample: The sample info get from table tb_task_sample
    :param d_analysis: The analysis dir
    :return res: The yaml config file for amplicon pipeline
    """
    # Prepare the info
    res = {}
    # 项目名称
    res["project"] = info_task[0]
    # 测序区域
    region_map = {"ITS": "ITS2", "18S": "V4"}
    info_type = info_task[20].strip().split(',')
    res["type"] = info_type[0]
    res["region"] = info_type[1] if info_type[0] == "16S" else region_map[info_type[0]]
    # 数据库
    db_map = {"16S": "greengenes", "18S": "silva", "ITS": "unite"}
    res["database"] = db_map[info_type[0]]
    # 样本路径
    d_data = os.path.join(d_analysis, "rawdata")
    os.makedirs(d_data, exist_ok=True)
    res["rawdata"] = d_data
    res["samples"] = linkdata(info_sample, d_data)
    # 数据类型
    res["data_type"] = "PE" if info_task[10].startswith("PE") else "SE"
    res["data_quality"] = 33
    # 分组信息
    info_group = get_group_info(info_sample)
    if info_group:
        res["groups"] = info_group
        res["group_order"] = list(info_group.keys())
    # 差异分析
    info_diff = get_diff_info(info_sample)
    if info_diff:
        res["group_diff"] = info_diff
    # 线程数，并行数
    res["threads"] = 8
    res["parallel"] = 8

    # 生成配置文件
    f_config = os.path.join(d_analysis, "config.yml")
    with open(f_config, 'w') as OUT:
        print(yaml.dump(res, sort_keys=False, allow_unicode=True), file=OUT)

    return f_config


####


def juage():
    db = pymysql.connect(
        host="172.16.0.18",
        user="vmtest",
        password="vmtest888",
        database="weiyuan",
        charset="utf8")

    # 使用 cursor() 方法创建一个游标对象 cursor
    cursor = db.cursor()

    # 分析中状态
    analyse = "ANALYSEING"
    sql_analyse = "select * from tb_analyse_task where status='%s'" % analyse

    cursor.execute(sql_analyse)
    arr_analyse = cursor.fetchone()

    # 排队状态
    queue = "NOT_START"
    # queue = "COMPLETE"
    sql_quese = "select * from tb_analyse_task where status='%s'" % queue

    cursor.execute(sql_quese)
    arr = cursor.fetchone()

    # 判断是否有任务在运行中，以及是否有任务在排队
    if not arr_analyse:
        print("没有任务处于运行中")
        if not arr:
            print("没有任务在排队")
        else:
            task_id = arr[0]  # 任务id，task_id

            # 样本信息表,提取task对应的样本id，以及样本的信息
            sql_task_sample = "select * from tb_task_sample where task_id='%s'" % task_id
            cursor.execute(sql_task_sample)
            cds = cursor.fetchall()
            df_all_sample = pd.DataFrame(list(cds))

            list1 = df_all_sample[7].tolist()
            print(list1)

            mkdir_shell = "mkdir -p /sdbb/Earth/%s /sdbb/Earth/Analysis/%s" % (task_id, task_id)
            os.system(mkdir_shell)

            sample = "/sdbb/Earth/%s/task_sample.csv" % task_id

            df_all_sample = df_all_sample[[7, 23, 24]]
            df_all_sample.to_csv(sample, header=False, index=False, sep='\t')
            # print (df_all_sample)

            log = "/sdbb/Earth/%s/log.txt" % task_id
            f = open(log, 'w')

            # 运行时，先把运行状态改成分析中
            sql_analyse = "update tb_analyse_task set status='ANALYSEING' where id='%s'" % task_id
            cursor.execute(sql_analyse)
            db.commit()

            type = arr[13]  # 类型
            sequencing_platform = arr[7]  # 测序平台
            created_library_strategy = arr[11]  # 建库策略
            pathogene = arr[12]  # 病原体

            if (type == "APPRAISAL_ANALYSIS"):
                if (sequencing_platform == "Nanopore"):
                    print("任务", task_id, time.asctime(), "开始nanopore鉴定流程")
                    shell = "perl /sdbb/share/pipeline/nanopore/bin/nano_identify.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" % (
                        task_id, sample, task_id)
                    print(shell)
                    os.system(shell)

                    return2mysql(list1=list1, cursor=cursor, task_id=task_id)

                    f.write("任务%s %s 完成nanopore鉴定流程" % (task_id, time.asctime()))
            elif (type == "NOVEL_CORONAVIRUS"):
                if (sequencing_platform == "Nanopore" and pathogene == "2019新型冠状病毒"):
                    print("任务", task_id, time.asctime(), "开始Nano新冠全基因组流程")
                    shell = "perl /sdbb/share/pipeline/nanopore/bin/nano_sars.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" % (
                        task_id, sample, task_id)
                    print(shell)
                    os.system(shell)

                    return2mysql(list1=list1, cursor=cursor, task_id=task_id)
                    f.write("任务%s %s 完成新冠全基因组流程" % (task_id, time.asctime()))
            elif (type == "PE150"):
                print("任务", task_id, time.asctime(), "开始宏基因PE150测序")
                # sh PE150.sh
                f.write("任务%s %s 完成宏基因PE150测序" % (task_id, time.asctime()))
            elif (type == "SIXTEEM_S"):
                print("任务", task_id, time.asctime(), "开始16S测序流程")
                # Add by renchaobo
                d_out = f"/sdbb/Earth/Analysis/{task_id}"
                f_config = config4amplicon(arr, cds, d_out)
                shell = f"python /sdbb/Earth/Develop/16S/main.py -c {f_config} -o {d_out} --run"
                os.system(shell)
                cmd = f"cp -rf {d_out}/{task_id} /sdbb/Earth/AnalysisResults/{task_id}\ncp {d_out}/{task_id}.zip /sdbb/Earth/AnalysisResults/{task_id}/{task_id}.zip"
                return2mysql(list1=list1, cursor=cursor, task_id=task_id)
                f.write("任务%s %s 完成16S测序流程" % (task_id, time.asctime()))
            elif (type == "DENOVO_JOINT"):
                if (created_library_strategy == "meta"):
                    print("任务", task_id, time.asctime(), "开始宏基因组组装流程")
                    shell = "perl /sdbb/bioinfor/lanlei/script/meta_denovo/bin/meta_denovo_and_anno.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" % (
                        task_id, sample, task_id)
                    print(shell)
                    os.system(shell)
                    # subprocess

                    f.write("任务'%s','%s'完成宏基因组组装流程\n" % (task_id, time.asctime()))
                elif (created_library_strategy == "bacteria"):
                    print("任务", task_id, time.asctime(), "开始宏基因组组装流程")
                    # sh DENOCO_bac.sh
                    f.write("任务'%s','%s'完成宏基因组组装流程\n" % (task_id, time.asctime()))
            else:
                # sh other.sh
                print("任务", task_id, time.asctime(), "未知流程")

            # result_dir = '/sdbb/bioinfor/EARTH/AnalysisResults/' + str(task_id)
            # sql_analyse = "update tb_task_sample set return_file_path='%s' where id='%s'" %(result_dir,task_id)
            # cursor.execute(sql_analyse)

            sql_end = "update tb_analyse_task set status='COMPLETE' where id='%s'" % task_id
            cursor.execute(sql_end)
            db.commit()

    else:
        print("有任务在运行中")

    cursor.close()  # 关闭游标
    db.close()  # 关闭数据库连接


def loop_func(func, second):
    # 每隔second秒执行func函数
    while True:
        func()
        time.sleep(second)


loop_func(juage, 10)
