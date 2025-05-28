import os
import time
import subprocess
import re
import logging

import pandas as pd
import pymysql
import yaml
from file.CheckWebInput import CheckWebInput
from file.get_16s_input import config4amplicon
from file.ei1010.common import dict2yaml  # dict转yaml
from file.ei1010 import pipe as mpipe


#############################
# mysql 操作
# db.cursor()
# execute() 执行
# fetchone() 方法获取单条数据
# fetchall() 方法获取多条数据
#############################

########################################################################################
########################################################################################
########################################################################################
# 返回，包括压缩包，html和fasta文件
def return2mysql(sample_lists,
                 cursor,
                 task_id,
                 zip_pattern=r"/sdbb/Earth/AnalysisResults/%s/%s.zip",
                 html_pattern=r"/file/weiyuan/AnalysisResults/%s/%s/index.html",
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
    for sample_id in sample_lists:
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
        cursor.execute(sql_url)
        # sequencing fasta
        if fa_pattern != "":
            # /sdbb/Earth/Analysis/2333/XG0120/4.consensus/XG0120.consensus.fa
            # /sdbb/Earth/Analysis/%s/%s/4.consensus/consensus.fa
            # fa_pattern = "/sdbb/Earth/AnalysisResults/%s/%s/2.assemble/%s.fasta"
            return_fasta = fa_pattern % (task_id, sample_id, sample_id)
            sql_fasta = "update tb_task_sample set return_fasta='%s' where task_id='%s' and sample_number='%s'" % (
                return_fasta, task_id, sample_id)
            cursor.execute(sql_fasta)

    chmod_shell = f'chmod -R 777 /sdbb/Earth/AnalysisResults/{task_id}'
    os.system(chmod_shell)


def return2mysql_16s(sample_lists,
                     cursor,
                     task_id,
                     zip_pattern=r"/sdbb/Earth/AnalysisResults/%s/%s.zip",
                     html_pattern=r"/file/weiyuan/AnalysisResults/%s/Upload/index.html",
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
    for sample_id in sample_lists:
        # zip
        return_file_path = zip_pattern % (task_id, task_id)
        sql_path = "update tb_task_sample set return_file_path='%s' where task_id='%s' and sample_number='%s'" % (
            return_file_path, task_id, sample_id)
        print(sql_path)
        cursor.execute(sql_path)
        # index html
        return_html_url = html_pattern % (task_id)
        sql_url = "update tb_task_sample set return_html_url='%s' where task_id='%s' and sample_number='%s'" % (
            return_html_url, task_id, sample_id)
        cursor.execute(sql_url)
        # sequencing fasta
        if fa_pattern != "":
            # /sdbb/Earth/Analysis/2333/XG0120/4.consensus/XG0120.consensus.fa
            # /sdbb/Earth/Analysis/%s/%s/4.consensus/consensus.fa
            return_fasta = fa_pattern % (task_id, sample_id)
            sql_fasta = "update tb_task_sample set return_fasta='%s' where task_id='%s' and sample_number='%s'" % (
                return_fasta, task_id, sample_id)
            cursor.execute(sql_fasta)
    chmod_shell = f'chmod -R 777 /sdbb/Earth/AnalysisResults/{task_id}'
    os.system(chmod_shell)


def return2mysql_typing(sample_lists,
                        cursor,
                        task_id,
                        zip_pattern=r"/sdbb/Earth/AnalysisResults/%s/%s.zip",
                        html_pattern=r"/file/weiyuan/AnalysisResults/%s/%s/index.html",):
    """
    输入: 1.样本名列表 [A1,A2,A3,...]
          2.MySQL游标对象
          3.任务ID
          3.单样本结果压缩包路径通配,第一个%s是任务ID,第二个%s是样本ID
          4.单样本网页结果路径通配 
          5.单样本结果FASTA序列路径,默认""无FASTA结果
    结果更新到MySQL数据库上,对应任务ID和样本ID的返回位置
    """
    for sample_id in sample_lists:
        # zip
        return_file_path = zip_pattern % (task_id, sample_id)
        sql_path = "update tb_task_sample set return_genotype_file_path='%s' where task_id='%s' and sample_number='%s'" % (
            return_file_path, task_id, sample_id)
        print(sql_path)
        cursor.execute(sql_path)
        # index html
        return_html_url = html_pattern % (task_id, sample_id)
        sql_url = "update tb_task_sample set return_genotype_html_url='%s' where task_id='%s' and sample_number='%s'" % (
            return_html_url, task_id, sample_id)
        cursor.execute(sql_url)

########################################################################################
########################################################################################
########################################################################################
# 任务管理。可以终止分析,并返回log文件等


def wrap_popen(cml: str, taskid, loghandle: str, cursor, timeout: int = 604800):
    """
    针对 Earth 对 subprocess.Popen 的包装
     参数 - cml:         生信流程的主脚本 shell command line
     参数 - taskid:      Earth 任务 ID
     参数 - loghandle:   记录文件 open 句柄
     参数 - cursor:      pymysql cursor 对象
     参数 - timeout:     超时设置, 单位秒
    """
    try:
        pp = subprocess.Popen(cml, shell=True, encoding="utf-8",
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        popen_perftime = time.perf_counter()
        while True:
            time.sleep(0.1)
            pp.poll()
            if pp.returncode == 0:
                break
            elif (pp.returncode is not None) and (pp.returncode != 0):
                break
            # 获取taskid对应的状态，CANCEL则终止进程
            tmp_sql = "select status from tb_analyse_task where task_id={}".format(
                taskid)
            cursor.execute(tmp_sql)
            current_status = cursor.fetchone()['status']
            if current_status == "CANCEL":
                pp.kill()
                break
            # 超过时间也会断掉
            if time.perf_counter() - popen_perftime > timeout:
                pp.kill()
                break
        # 运行记录 & Popen 返回码
        if pp.returncode is None:
            loghandle.write("终止运行!\n")
        else:
            loghandle.write("[STDOUT]\n" + pp.stdout.read())
            loghandle.write("\n[STDERR]\n" + pp.stderr.read())
        return pp.returncode
    except Exception as e:
        print(e)


def wrap_popen_no_stderr(cml: str, taskid, loghandle: str, cursor, timeout: int = 604800):
    """
    针对 Earth 对 subprocess.Popen 的包装
    # rcb
    # 消除死锁

     参数 - cml:         生信流程的主脚本 shell command line
     参数 - taskid:      Earth 任务 ID
     参数 - loghandle:   记录文件 open 句柄
     参数 - cursor:      pymysql cursor 对象
     参数 - timeout:     超时设置, 单位秒
    """
    try:
        pp = subprocess.Popen(cml, shell=True, encoding="utf-8",
                              stdout=loghandle, stderr=loghandle)
        popen_perftime = time.perf_counter()
        while True:
            time.sleep(0.1)
            pp.poll()
            if pp.returncode == 0:
                break
            elif (pp.returncode is not None) and (pp.returncode != 0):
                break
            # 获取taskid对应的状态，CANCEL则终止进程
            tmp_sql = "select status from tb_analyse_task where task_id={}".format(
                taskid)
            cursor.execute(tmp_sql)
            current_status = cursor.fetchone()['status']
            if current_status == "CANCEL":
                pp.kill()
                break
            # 超过时间也会断掉
            if time.perf_counter() - popen_perftime > timeout:
                pp.kill()
                break
        # 运行记录 & Popen 返回码
        if pp.returncode is None:
            loghandle.write("终止运行!\n")
        else:
            loghandle.write("[STDOUT]\n" + pp.stdout.read())
            loghandle.write("\n[STDERR]\n" + pp.stderr.read())
        return pp.returncode
    except Exception as e:
        print(e)
########################################################################################
########################################################################################
########################################################################################
# 返回log文件，不同的工作流返回的不一样，有的返回zip，有的直接返回log文件


def log_return(workflow_name, task_id):
    if (workflow_name == "16S测序" or workflow_name == "Blast"):
        cmd = "cp -rf /sdbb/Earth/Analysis/%s/log.txt /sdbb/Earth/AnalysisResults/%s/" % (
            task_id, task_id)
        os.system(cmd)
        return "/sdbb/Earth/AnalysisResults/%s/log.txt" % task_id
    else:
        cmd = "cp -rf /sdbb/Earth/Analysis/%s/%s_log.zip /sdbb/Earth/AnalysisResults/%s/" % (
            task_id, task_id, task_id)
        os.system(cmd)
        return "/sdbb/Earth/AnalysisResults/%s/%s_log.zip" % (task_id, task_id)

########################################################################################
########################################################################################
########################################################################################
# 220330 mengxf 没生成 index.html 判断为运行失败 ERROR


def check_results(sample_lists,
                  task_id,
                  workflow_name):
    """
    输入: 1.样本名列表 [A1,A2,A3,...]
          2.任务ID
          3.工作流名
    返回: 返回与"wpopen"一致的返回码, 成功为0, 失败为 1
    """
    # 一个index.html
    if workflow_name in ["16S测序", "宏基因组组间分析", "泛基因组分析"]:
        # /sdbb/Earth/AnalysisResults/1683/Upload/index.html
        index_html = f"/sdbb/Earth/AnalysisResults/{task_id}/Upload/index.html"
        if not os.path.isfile(index_html):
            return 1
    elif workflow_name == "溯源进化树分析":
        index_html = f"/sdbb/Earth/AnalysisResults/{task_id}/index.html"
        if not os.path.isfile(index_html):
            return 1
    # 每个index.html
    else:
        for sample in sample_lists:
            # /sdbb/Earth/AnalysisResults/1438/s1/index.html
            index_html = f"/sdbb/Earth/AnalysisResults/{task_id}/{sample}/index.html"
            if not os.path.isfile(index_html):
                return 1
    return 0

########################################################################################
########################################################################################
########################################################################################


def bwa_index(ref_index, ref_seq, ref_genome_id, cursor):
    base_name = os.path.basename(ref_seq)
    shell = f"cp -R {ref_seq} /sdbb/Earth/database/ref_genome/\n"
    if (ref_index == 0):
        if ('gz' in ref_seq):
            shell += f"gzip -d /sdbb/Earth/database/ref_genome/{base_name}\n"
            # 1.fa.gz 1.fa
            base_name = base_name.rstrip(".gz")
        # else:
        #     pass
        shell += f"bwa index /sdbb/Earth/database/ref_genome/{base_name}"
        print(shell)
        try:
            os.system(shell)
            sql_index = "update tb_dict_ref_gene set ref_index='1' where ref_id='%s'" % ref_genome_id
            cursor.execute(sql_index)
        except Exception as e:
            print(e)

    return "/sdbb/Earth/database/ref_genome/%s" % base_name


########################################################################################
########################################################################################
########################################################################################

def juage():
    db = pymysql.connect(
        # host="172.16.0.18",
        host="localhost",
        user="vmtest",
        password="vmtest888",
        database="weiyuan",
        autocommit=True,
        charset="utf8")

    # db = pymysql.connect(
    # host="testing.idaoben.com",
    # port=3524,
    # user="weiyuan",
    # password="0BUsk5ZWvf7y3KYt",
    # database="weiyuan",
    # autocommit=True,
    # charset="utf8")

    # 使用 cursor() 方法创建一个游标对象 cursor
    cursor = db.cursor(cursor=pymysql.cursors.DictCursor)

    # 分析中状态
    analyse = "ANALYSEING"
    sql_analyse = "select * from tb_analyse_task where status='%s'" % analyse

    cursor.execute(sql_analyse)
    arr_analyse = cursor.fetchone()

    # 排队状态
    queue = "NOT_START"
    sql_quese = "select * from tb_analyse_task where status='%s'" % queue

    cursor.execute(sql_quese)
    arr = cursor.fetchone()

    # 判断是否有任务在运行中，以及是否有任务在排队
    if not arr_analyse:
        print("没有任务处于运行中")
        if not arr:
            print("没有任务在排队中")
        else:
            task_id = arr['task_id']
            print("任务ID:" + str(task_id))
            # 样本信息表,提取task对应的样本id，以及样本的信息
            sql_task_sample = "select * from tb_task_sample where task_id='%s'" % task_id
            cursor.execute(sql_task_sample)
            cds = cursor.fetchall()
            df_all_sample = pd.DataFrame(list(cds))

            # 第7列为样本id，存在sample_lists上面，后续根据列表来返回html等。
            sample_lists = df_all_sample['sample_number'].tolist()  # 样本名
            print(sample_lists)

            # 三个目录，分别存放原始数据，分析过程，分析结果
            rm_shell = f"rm /sdbb/Earth/Analysis/{task_id} -rf"
            os.system(rm_shell)
            mkdir_shell = "mkdir -p /sdbb/Earth/Rawdata/%s /sdbb/Earth/Analysis/%s /sdbb/Earth/AnalysisResults/%s" % (
                task_id, task_id, task_id)
            os.system(mkdir_shell)
            chmod_shell = "chmod -R 777 /sdbb/Earth/Rawdata/%s  /sdbb/Earth/AnalysisResults/%s" % (
                task_id, task_id)
            os.system(chmod_shell)

            # 第7,23,24列分别为name、fq1、fq2，大部分流程只需要这三个数据
            sample = "/sdbb/Earth/Analysis/%s/task_sample.csv" % (task_id)
            df_all_sample_filter = df_all_sample[[
                'sample_number', 'sequencing_data1', 'sequencing_data2']]
            df_all_sample_filter.to_csv(
                sample, header=False, index=False, sep='\t')

            # conf_meta = "/sdbb/Earth/Analysis/%s/group.csv" %(task_id)
            # df_meta_mapping = df_all_sample[['sample_number','sequencing_data1','sequencing_data2','label','label_compare']]
            # df_meta_mapping.to_csv(conf_meta,header=False,index=False,sep='\t')

            # 存放log文件
            log = "/sdbb/Earth/AnalysisResults/%s/log.txt" % task_id
            f = open(log, 'w')

            # 运行时，先把运行状态改成分析中
            sql_analyse = "update tb_analyse_task set status='ANALYSEING' where task_id='%s'" % task_id
            cursor.execute(sql_analyse)
            db.commit()

            sars_database = arr["sars_cov_2_database"]  # 是否跑新冠全球库
            ref_genome_id = arr["reference_genome"]  # 参考基因组ID
            tree_data_id = arr["tree_data"]  # 进化树ID
            host_id = arr["host_id"]  # 宿主ID
            genome_typing_database_id = arr["genome_typing_database"]  # 分型库id
            mars_result_record = arr["mars_result_record"]  # mars记录的id
            # mars源或者外部导入
            new_pathogen_data_source = arr['new_pathogen_data_source']
            qc_parameter = arr['qc_parameter']  # QC流程excel文件
            pe_or_se = arr["pe_or_se"]
            sequence = ''  # 单端或者双端
            if (pe_or_se == "单端"):
                sequence = 'SE'
            elif (pe_or_se == "双端"):
                sequence = 'PE'

            # 参考基因组路径
            if (ref_genome_id != ""):
                sql_ref = "select * from tb_dict_ref_gene where ref_id='%s'" % ref_genome_id
                cursor.execute(sql_ref)
                cds_ref = cursor.fetchone()
                ref_seq = cds_ref['ref_seq']
                ref_index = cds_ref['ref_index']
                new_ref_path = bwa_index(
                    ref_index, ref_seq, ref_genome_id, cursor)
                # print (ref_seq)
                print("参考基因组位置:" + new_ref_path)

            # 进化树路径
            if (tree_data_id != ""):
                sql_tree = "select * from tb_dict_evolutionary_tree where tree_id='%s'" % tree_data_id
                cursor.execute(sql_tree)
                cds_tree = cursor.fetchone()
                tree_ref_fa = cds_tree['tree_ref_fa']
                print("进化树fa路径:" + tree_ref_fa)

            # 宿主名,中文
            if (host_id):
                sql_host = "select * from tb_dict_host where host_id='%s'" % host_id
                cursor.execute(sql_host)
                cds_host = cursor.fetchone()
                host_name = cds_host['host_name']
                print("宿主:" + host_name)

            # 分型库名字，中文
            if (genome_typing_database_id):
                sql_typing = "select * from tb_dict_genome_type where pathogen_id='%s'" % genome_typing_database_id
                cursor.execute(sql_typing)
                cds_typing = cursor.fetchone()
                pathogen_name = cds_typing['pathogen_name']
                print("分型库病原名:" + pathogen_name)

            # mars记录
            if (mars_result_record != ""):
                mars_library = ''
                mars_sublibrary = ''
                mars_pathogen = ''
                if (re.match(',', mars_result_record)):
                    sql_typing = "select * from tb_dict_mars where mars_id='%s'" % mars_result_record
                    cursor.execute(sql_typing)
                    cds_mars = cursor.fetchone()
                    mars_library = cds_mars['mars_library']
                    mars_sublibrary = cds_mars['mars_sublibrary']
                    mars_pathogen = cds_mars['pathogen_scientific_name']
                    print("mars文库:" + mars_library)
                    print("mars子文库:" + mars_sublibrary)
                    print("mars病原体:" + mars_pathogen)
                else:
                    print(f'mars_result_record{mars_result_record}')
                    all_ids = mars_result_record.split(',')
                    pathogens = []
                    for mars_id in all_ids:
                        print(f'mars_id{mars_id}')
                        sql_typing = "select * from tb_dict_mars where mars_id='%s'" % mars_id
                        cursor.execute(sql_typing)
                        cds_mars = cursor.fetchone()
                        mars_library = cds_mars['mars_library']
                        mars_sublibrary = cds_mars['mars_sublibrary']
                        pathogen_ma = cds_mars['pathogen_scientific_name']
                        print(pathogen_ma)
                        if (pathogen_ma):
                            pathogens.append(pathogen_ma)
                        else:
                            pathogens.append('unmap')
                        print("mars文库:" + mars_library)
                        print("mars子文库:" + mars_sublibrary)
                    print(pathogens)
                    mars_pathogen = ','.join(pathogens)
                    print("mars病原体:" + mars_pathogen)

            # 执行检查脚本
            check_web_input = CheckWebInput(task_id=task_id, cursor=cursor)
            dict2yaml(check_web_input.indict,
                      f"/sdbb/Earth/Analysis/{task_id}/task_samples.yaml")

            # 返回error—message，直接返回的就是字符串
            error_message = check_web_input.error_message
            # 如果error_message为空，正常执行后续流程，如果不为空，则返回ERROE状态，并将error_message写入到tb_analyse_task表中的error_message
            if (error_message != ""):
                err = "%s-错误信息为:%s" % (task_id, error_message)

                sql_error = "update tb_analyse_task set status='ERROR' where task_id='%s'" % task_id
                cursor.execute(sql_error)
                db.commit()

                sql_error = "update tb_analyse_task set error_message='%s' where task_id='%s'" % (
                    error_message[:255], task_id)
                try:
                    cursor.execute(sql_error)
                except:
                    pass

                db.commit()

                return None
            ###############################################################################################################
            ###############################################################################################################
            else:  # 前面检查没问题，走正常工作流
                # check_web_input返回工作流
                workflow = check_web_input.workflow
                dict2yaml(check_web_input.indict,
                          f"/sdbb/Earth/Analysis/{task_id}/task_samples.yaml")
                print("工作流程为:" + workflow)
                
                popen_returncode = None
                # 判断进入的工作流
                ###############################################################################################################
                if (workflow == "nanopore鉴定"):
                    shell = "perl /sdbb/share/pipeline/nanopore/bin/nano_identify.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" % (
                        task_id, sample, task_id)
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    print(popen_returncode)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "kraken鉴定"):
                    shell = f"perl /sdbb/share/pipeline/Species_identification/bin/species_identification.pl -fq_input {sample} -sequence {sequence} -workdir /sdbb/Earth/Analysis/{task_id} -task_id {task_id} -run"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    print(popen_returncode)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "新发病原分析"):
                    if (new_pathogen_data_source == "Mars记录"):
                        input_shell = f"perl  /sdbb/share/pipeline/new_pathogen/bin/get_mars_input_fq.pl -library {mars_library} -sublibrary {mars_sublibrary} -workdir /sdbb/Earth/Analysis/{task_id} "
                        os.system(input_shell)
                        print(input_shell)

                        shell = f"perl /sdbb/share/pipeline/new_pathogen/bin/new_pathogen.pl -mars_or_out mars -library {mars_library} -sublibrary {mars_sublibrary} -workdir /sdbb/Earth/Analysis/{task_id} -task_id {task_id}"
                    else:
                        shell = f"perl /sdbb/share/pipeline/new_pathogen/bin/new_pathogen.pl -mars_or_out out -fq_input {sample} -workdir /sdbb/Earth/Analysis/{task_id} -task_id {task_id}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                elif (workflow == "mars导入全基因组"):
                    # mars_library = cds_mars['mars_library']
                    # mars_sublibrary = cds_mars['mars_sublibrary']
                    # mars_pathogen = cds_mars['pathogen_scientific_name']

                    shell = f"perl /sdbb/share/pipeline/mars_wgs/bin/mars_wgs.pl -library {mars_library} -sublibrary {mars_sublibrary} -pathogen_name {mars_pathogen} -mars_oudir /sdbb/Data_back/IDseqV2/ -workdir /sdbb/Earth/Analysis/{task_id} -task_id {task_id}\n"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "Nano新冠全基因组"):
                    if (sars_database == "否"):
                        shell = "perl /sdbb/share/pipeline/nanopore/bin/nano_sars.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s -database n" % (
                            task_id, sample, task_id)
                    elif (sars_database == "是"):
                        shell = "perl /sdbb/share/pipeline/nanopore/bin/nano_sars.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s -database y" % (
                            task_id, sample, task_id)
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "二代测序新冠扩增子分析"):
                    # 240228 mengxf
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.ngs_sars2_wga()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "二代测序流感扩增子分析"):
                    # 20220627 RCB
                    cmd = f"python /sdbb/share/pipeline/FluAmplicon/main.py --venus /sdbb/Earth/Analysis/{task_id}/task_samples.yaml"
                    popen_returncode = wrap_popen(cmd, task_id, f, cursor)
                    if popen_returncode == 0:
                        cmd = f"cd /sdbb/Earth/Analysis/{task_id};zip -r {task_id}_log.zip Shell >/dev/null 2>&1;cp {task_id}_log.zip -rf /sdbb/Earth/AnalysisResults/{task_id};cp -rf Upload/* /sdbb/Earth/AnalysisResults/{task_id}"
                        subprocess.run(cmd, shell=True)
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                    # shell = "bash /sdbb/share/pipeline/Sars_Cov2_Amplicon/sars_cov2_wrap_for_earthFQ.sh %s %s" %(sample,task_id)

                    # popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    # if popen_returncode == 0:
                        # return2mysql在fa_pattern不为空的时候，会自动补齐几个参数
                        # return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id,fa_pattern="/sdbb/Earth/AnalysisResults/%s/%s/4.consensus/%s.consensus.fa")
                ###############################################################################################################
                elif (workflow == "IDseqComplete_超快多重_PCR_新冠全基因组建库试剂盒"):
                    # 221102 mengxf 公司新冠试剂盒专用流程
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.ngs_sars2_wga_idseqcomplete()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "超多重PCR扩增病原鉴定"):
                    cmd = f"python /sdbb/share/pipeline/ampliconAnalysis/main.py -i /sdbb/Earth/Analysis/{task_id}/task_samples.yaml  -o /sdbb/Earth/Analysis/{task_id} -u /sdbb/Earth/AnalysisResults/{task_id}"
                    print(cmd)
                    popen_returncode = wrap_popen(cmd, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "新冠-污水监测"):
                    # 230417 mengxf ILMN新冠污水监测流程
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.ngs_sars_wastewater_surveillance()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "nanopore新冠-污水监测"):
                    # 230425 mengxf ONT新冠污水监测流程
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.ont_sars_wastewater_surveillance()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "二代测序其他扩增子分析"):
                    # 240228 mengxf
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.ngs_amplicon_general()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "单病毒三代denovo组装"):
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    # 无参
                    if not check_web_input.indict["reference_genome"]:
                        mypipe.ont_virus_assembly()
                    # 有参
                    else:
                        mypipe.ont_virus_map()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "单细菌三代denovo组装"):
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.ont_bacteria_assembly()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "单病毒混合denovo组装"):
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.hyb_virus_assembly()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "单细菌混合denovo组装"):
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.hyb_bacteria_assembly()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "单病毒的有参组装"):
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.ont_virus_map()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "单细菌denovo组装"):
                    if not tree_data_id:
                        tree_ref_fa = "blank"
                    if not genome_typing_database_id:
                        pathogen_name = "blank"

                    ref_seq = ref_seq if ref_genome_id else 'None'
                    shell = f"python /sdbb/share/pipeline/single_denovo/parallel.py -kind bacteria -table /sdbb/Earth/Analysis/{task_id}/task_samples.yaml -out /sdbb/Earth/Analysis/{task_id} -upload /sdbb/Earth/AnalysisResults/{task_id} -tree {tree_ref_fa}  -typing {pathogen_name} -ref {ref_seq}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists, cursor=cursor, task_id=task_id,
                                     fa_pattern="/sdbb/Earth/AnalysisResults/%s/%s/2.ass/%s_scaffolds.fasta")
                        return2mysql_typing(sample_lists=sample_lists, cursor=cursor, task_id=task_id, zip_pattern="/sdbb/Earth/AnalysisResults/%s/%s/Tree_Typing.zip",
                                            html_pattern="/file/weiyuan/AnalysisResults/%s/%s/Tree_Typing/index.html")
                ###############################################################################################################
                elif (workflow == "单真菌denovo组装"):
                    if not tree_data_id:
                        tree_ref_fa = "blank"
                    if not genome_typing_database_id:
                        pathogen_name = "blank"
                    ref_seq = ref_seq if ref_genome_id else 'None'
                    shell = f"python /sdbb/share/pipeline/single_denovo/parallel.py -kind fungi -table /sdbb/Earth/Analysis/{task_id}/task_samples.yaml -out /sdbb/Earth/Analysis/{task_id} -upload /sdbb/Earth/AnalysisResults/{task_id} -tree {tree_ref_fa}   -typing {pathogen_name} -ref {ref_seq}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists, cursor=cursor, task_id=task_id,
                                     fa_pattern="/sdbb/Earth/AnalysisResults/%s/%s/2.ass/%s_scaffolds.fasta")
                        return2mysql_typing(sample_lists=sample_lists, cursor=cursor, task_id=task_id, zip_pattern="/sdbb/Earth/AnalysisResults/%s/%s/Tree_Typing.zip",
                                            html_pattern="/file/weiyuan/AnalysisResults/%s/%s/Tree_Typing/index.html")
                ###############################################################################################################
                elif (workflow == "单病毒denovo组装"):
                    if not tree_data_id:
                        tree_ref_fa = "blank"
                    if not genome_typing_database_id:
                        pathogen_name = "blank"
                    ref_seq = ref_seq if ref_genome_id else 'None'
                    shell = f"python /sdbb/share/pipeline/single_denovo/parallel.py -kind virus -table /sdbb/Earth/Analysis/{task_id}/task_samples.yaml -out /sdbb/Earth/Analysis/{task_id} -upload /sdbb/Earth/AnalysisResults/{task_id} -tree {tree_ref_fa}  -typing {pathogen_name}  -ref {ref_seq}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists, cursor=cursor, task_id=task_id,
                                     fa_pattern="/sdbb/Earth/AnalysisResults/%s/%s/2.align/%s_consensus.fa")
                        return2mysql_typing(sample_lists=sample_lists, cursor=cursor, task_id=task_id, zip_pattern="/sdbb/Earth/AnalysisResults/%s/%s/Tree_Typing.zip",
                                            html_pattern="/file/weiyuan/AnalysisResults/%s/%s/Tree_Typing/index.html")
                ###############################################################################################################
                elif (workflow == "新冠宏基因组PE150拼接"):
                    shell = f"python /sdbb/share/pipeline/nCov2019_meta_assemble/parallel.py -slist {sample} -out /sdbb/Earth/Analysis/{task_id} -upload /sdbb/Earth/AnalysisResults/{task_id}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "新冠宏基因组SE50拼接"):
                    shell = f"python /sdbb/share/pipeline/nCov2019_meta_assemble/parallel.py -slist {sample} -out /sdbb/Earth/Analysis/{task_id} -upload /sdbb/Earth/AnalysisResults/{task_id}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "其他宏基因组PE150拼接"):
                    shell = f"perl /sdbb/share/pipeline/ncov2019_blast_denovo/bin/ncov2019_blast_denovo.pl -workdir /sdbb/Earth/Analysis/%s -sequence PE -fq_input %s -task_id %s -ref %s" % (
                        task_id, sample, task_id, new_ref_path)

                    print(shell)
                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "宏基因组无参组装"):
                    shell = "perl /sdbb/share/pipeline/meta_denovo/bin/meta_denovo_and_anno.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" % (
                        task_id, sample, task_id)
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "mars导入全基因组"):
                    mars_library = cds_mars['mars_library']
                    mars_sublibrary = cds_mars['mars_sublibrary']
                    mars_pathogen = cds_mars['pathogen_scientific_name']

                    shell = f"perl /sdbb/share/pipeline/mars_wgs/bin/mars_wgs.pl -library {mars_library} -sublibrary {mars_sublibrary} -pathogen_name {mars_pathogen} -mars_oudir /sdbb/Data_back/IDseqV2/ -workdir /sdbb/Earth/Analysis/{task_id} -task_id {task_id}\n"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "病原分型分析"):
                    shell = f"python /sdbb/share/pipeline/TypeTracePlugin/main.py --config /sdbb/Earth/Analysis/{task_id}/task_samples.yaml >/sdbb/Earth/AnalysisResults/{task_id}/log1.txt 2>&1"
                    print(shell)

                    # popen_returncode = wrap_popen_no_stderr(shell, task_id, f, cursor)
                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql_typing(
                            sample_lists=sample_lists, cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "耐药毒力分析"):
                    d_analy = f"/sdbb/Earth/Analysis/{task_id}"
                    d_result = f"/sdbb/Earth/AnalysisResults/{task_id}"
                    shell = f"python3 /sdbb/share/pipeline/Virus_drug_Plugin/multi_main.py -table /sdbb/Earth/Analysis/{task_id}/task_samples.yaml -out {d_analy} -upload {d_result}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "溯源进化树分析"):
                    mypipe = mpipe.EI1010(check_web_input.indict, cursor)
                    mypipe.trace_tree()
                    popen_returncode = mypipe.popen_returncode
                ###############################################################################################################
                elif (workflow == "宏基因组组间分析"):
                    conf_meta = "/sdbb/Earth/Analysis/%s/raw_input.csv" % (
                        task_id)
                    df_meta_mapping = df_all_sample[[
                        'sample_number', 'sequencing_data1', 'sequencing_data2', 'label', 'label_compare']]
                    df_meta_mapping.to_csv(
                        conf_meta, header=False, index=False, sep='\t')
                    shell_config = f"perl /sdbb/share/pipeline/meta_bin/bin/deal_mysql_csv2config.pl {conf_meta} /sdbb/Earth/Analysis/{task_id}"
                    try:
                        os.system(shell_config)
                    except Exception as e:
                        print(e)
                    shell = f"perl /sdbb/share/pipeline/meta_bin/bin/meta_mapping.pl -fq_input {sample} -group /sdbb/Earth/Analysis/{task_id}/group.list -color /sdbb/Earth/Analysis/{task_id}/color.list -diff /sdbb/Earth/Analysis/{task_id}/diff.list -workdir /sdbb/Earth/Analysis/{task_id}/ -task_id {task_id}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql_16s(
                            sample_lists=sample_lists, cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "16S测序"):
                    d_out = f"/sdbb/Earth/Analysis/{task_id}"
                    f_config = config4amplicon(arr, cds, d_out)
                    shell = f"python /sdbb/share/pipeline/16S/main.py -c {f_config} -o {d_out} --run >/sdbb/Earth/AnalysisResults/{task_id}/log1.txt 2>&1"

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql_16s(
                            sample_lists=sample_lists, cursor=cursor, task_id=task_id)

                    cmd = f"cp -rf {d_out}/Upload /sdbb/Earth/AnalysisResults/{task_id}\ncp {d_out}/{task_id}.zip /sdbb/Earth/AnalysisResults/{task_id}/{task_id}.zip"
                    os.system(cmd)

                ###############################################################################################################
                elif (workflow == "fastq质控"):
                    shell = f"""python /sdbb/share/pipeline/QC/bin/get_yaml_from_xlsx.py {qc_parameter} /sdbb/Earth/Analysis/{task_id}/qc.yaml
python /sdbb/share/pipeline/QC/bin/qc_fastp.py --input {sample} --f_yaml /sdbb/Earth/Analysis/{task_id}/qc.yaml --outdir /sdbb/Earth/Analysis/{task_id} --sequence {sequence} --run --task_id {task_id}
"""
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "泛基因组分析"):
                    shell = f"""
python /sdbb/share/pipeline/core_gene/bin/get_core_Input.py --f_yaml /sdbb/Earth/Analysis/{task_id}/task_samples.yaml --output /sdbb/Earth/Analysis/{task_id}/core_fa.csv
python /sdbb/share/pipeline/core_gene/bin/core_gene.py --input /sdbb/Earth/Analysis/{task_id}/core_fa.csv --outdir /sdbb/Earth/Analysis/{task_id} --run --task_id {task_id}
"""
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql_16s(
                            sample_lists=sample_lists, cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "蛋白质组学分析"):
                    shell = f"python3 /sdbb/share/pipeline/ProteomicsPlugin/main.py --solar_yaml /sdbb/Earth/Analysis/{task_id}/task_samples.yaml --out /sdbb/Earth/Analysis/{task_id} --upload /sdbb/Earth/AnalysisResults/{task_id}"
                    print(shell)
                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "代谢组学分析"):
                    pass
                ###############################################################################################################
                elif (workflow == "Blast"):
                    d_analy = f"/sdbb/Earth/Analysis/{task_id}"
                    d_result = f"/sdbb/Earth/AnalysisResults/{task_id}"
                    shell = f"python3 /sdbb/share/pipeline/BlastPlugin/multi_blast.py -data /sdbb/Earth/Analysis/{task_id}/task_samples.yaml -out {d_analy} -upload {d_result}"

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                    cmd = f"cp -rf {d_analy}/Upload {d_result};cp -rf {d_analy}/Upload.zip {d_result}/{task_id}.zip"
                    os.system(cmd)

                ###############################################################################################################
                elif (workflow == "去宿主"):
                    d_analy = f"/sdbb/Earth/Analysis/{task_id}"
                    d_result = f"/sdbb/Earth/AnalysisResults/{task_id}"
                    shell = f"python3 /sdbb/share/pipeline/Host_Plugin/multi_trim_host.py -table {sample} -host {host_name} -out {d_analy} -upload {d_result}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)

                ###############################################################################################################
                elif (workflow == "checkM工具"):
                    d_analy = f"/sdbb/Earth/Analysis/{task_id}"
                    d_result = f"/sdbb/Earth/AnalysisResults/{task_id}"
                    shell = f"python3 /sdbb/share/pipeline/CheckmPlugin/main.py -y /sdbb/Earth/Analysis/{task_id}/task_samples.yaml --upload {d_result}"
                    print(shell)

                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql(sample_lists=sample_lists,
                                     cursor=cursor, task_id=task_id)
                ###############################################################################################################
                elif (workflow == "数据拆分"):
                    pass

                ###############################################################################################################
                sql_end = "update tb_analyse_task set status='COMPLETE' where task_id='%s'" % task_id
                cursor.execute(sql_end)
                db.commit()

                # 检查结果是否有生成html，如果一个task—_id里面有至少一个样本没有html，check_result == 1，报错
                check_result = check_results(
                    sample_lists=sample_lists, task_id=task_id, workflow_name=workflow)
                # 最后
                if popen_returncode is None:
                    # 终止
                    f.write("任务%s %s 终止%s流程" %
                            (task_id, time.asctime(), workflow))
                    sql_end = "update tb_analyse_task set status='CANCEL' where task_id='%s'" % task_id
                elif (popen_returncode == 0 and check_result == 0):  # 全部样本都有生成html
                    # 完成
                    sql_end = "update tb_analyse_task set status='COMPLETE' where task_id='%s'" % task_id
                    f.write("任务%s %s 完成%s流程" %
                            (task_id, time.asctime(), workflow))
                else:  # popen接收到的错误，以及不是所有样本都有生成html
                    error_message_return = log_return(workflow, task_id)
                    sql_err = "update tb_analyse_task set error_message='%s' where task_id='%s'" % (
                        error_message_return, task_id)
                    cursor.execute(sql_err)
                    print(f"popen_returncode: {popen_returncode}")
                    print("sql_err" + sql_err)
                    db.commit()

                    # 异常中断
                    sql_end = "update tb_analyse_task set status='ERROR' where task_id='%s'" % task_id
                    f.write("任务%s %s 异常中断%s流程" %
                            (task_id, time.asctime(), workflow))

                print("sql_end:" + sql_end)
                cursor.execute(sql_end)
                db.commit()
###############################################################################################################
    else:
        print("有任务在运行中")
        pass

    cursor.close()  # 关闭游标
    db.close()  # 关闭数据库连接


def loop_func(func, second):
 # 每隔second秒执行func函数
    while True:
        func()
        time.sleep(second)


loop_func(juage, 10)
