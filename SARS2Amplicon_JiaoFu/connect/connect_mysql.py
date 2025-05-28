from email import message
import os
from sqlite3 import connect
import time
import subprocess
# from subprocess import run

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

########################################################################################
########################################################################################
########################################################################################
# 返回，包括压缩包，html和fasta文件
def return2mysql(sample_lists,
                cursor,
                task_id,
                zip_pattern = r"/sdbb/Earth/AnalysisResults/%s/%s.zip", 
                html_pattern = r"http://172.16.0.18:8080/file/weiyuan/AnalysisResults/%s/%s/index.html",
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
        return_file_path = zip_pattern %(task_id,sample_id)
        sql_path = "update tb_task_sample set return_file_path='%s' where task_id='%s' and sample_number='%s'" %(return_file_path,task_id,sample_id)
        print (sql_path)
        cursor.execute(sql_path)
        # index html
        return_html_url = html_pattern %(task_id,sample_id)
        sql_url = "update tb_task_sample set return_html_url='%s' where task_id='%s' and sample_number='%s'" %(return_html_url,task_id,sample_id)
        print (sql_url)
        cursor.execute(sql_url)
        # sequencing fasta
        if fa_pattern != "":
            # /sdbb/Earth/Analysis/2333/XG0120/4.consensus/XG0120.consensus.fa
            # /sdbb/Earth/Analysis/%s/%s/4.consensus/consensus.fa
            return_fasta = fa_pattern %(task_id,sample_id,sample_id)
            sql_fasta = "update tb_task_sample set return_fasta='%s' where task_id='%s' and sample_number='%s'" %(return_fasta,task_id,sample_id)
            cursor.execute(sql_fasta)


def return2mysql_16s(sample_lists,
                    cursor,
                    task_id,
                    zip_pattern = r"/sdbb/Earth/AnalysisResults/%s/%s.zip", 
                    html_pattern = r"/file/weiyuan/AnalysisResults/%s/Upload/index.html",
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
        return_file_path = zip_pattern %(task_id,task_id)
        sql_path = "update tb_task_sample set return_file_path='%s' where task_id='%s' and sample_number='%s'" %(return_file_path,task_id,sample_id)
        print (sql_path)
        cursor.execute(sql_path)
        # index html
        return_html_url = html_pattern %(task_id)
        sql_url = "update tb_task_sample set return_html_url='%s' where task_id='%s' and sample_number='%s'" %(return_html_url,task_id,sample_id)
        print (sql_url)
        cursor.execute(sql_url)
        # sequencing fasta
        if fa_pattern != "":
            # /sdbb/Earth/Analysis/2333/XG0120/4.consensus/XG0120.consensus.fa
            # /sdbb/Earth/Analysis/%s/%s/4.consensus/consensus.fa
            return_fasta = fa_pattern %(task_id,sample_id)
            sql_fasta = "update tb_task_sample set return_fasta='%s' where task_id='%s' and sample_number='%s'" %(return_fasta,task_id,sample_id)
            cursor.execute(sql_fasta)


########################################################################################
########################################################################################
########################################################################################

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
            subprocess.run(cmd, shell=True)
        elif sample_type == "SE50":
            cmd = f"ln -sf {record[23]} {d_rawdata}/{sample_name}_1.fq.gz"
            subprocess.run(cmd, shell=True)
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

########################################################################################
########################################################################################
########################################################################################

#[20220311 mengxf]######################################################################

# 任务管理。可以终止分析,并返回log文件等
def wrap_popen(cml:str, taskid, loghandle:str, cursor, timeout:int=604800):
    """
    针对 Earth 对 subprocess.Popen 的包装
     参数 - cml:         生信流程的主脚本 shell command line
     参数 - taskid:      Earth 任务 ID
     参数 - loghandle:   记录文件 open 句柄
     参数 - cursor:      pymysql cursor 对象
     参数 - timeout:     超时设置, 单位秒
    """
    pp = subprocess.Popen(cml, shell=True, encoding="utf-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print("Popen - 后台运行中")
    popen_perftime = time.perf_counter()
    while True:
        time.sleep(0.1)
        pp.poll()
        if pp.returncode == 0:
            print("Popen - 顺利执行！")
            break
        elif (pp.returncode is not None) and (pp.returncode != 0):
            print("Popen - 错误执行！")
            break
        # 获取taskid对应的状态，CANCEL则终止进程
        tmp_sql = "select status from tb_analyse_task where id={}".format(taskid)
        cursor.execute(tmp_sql)
        current_status = cursor.fetchone()[0]
        if current_status == "CANCEL":
            pp.kill()
            print("Popen - 接到终止信号")
            break
        # 超过时间也会断掉
        if time.perf_counter() - popen_perftime > timeout:
            pp.kill()
            print("Popen - 超时.设定: {}秒".format(timeout))
            break
    # 运行记录 & Popen 返回码
    if pp.returncode is None:
        loghandle.write("终止运行!\n")
    else:
        loghandle.write("[STDOUT]\n" + pp.stdout.read())
        loghandle.write("\n[STDERR]\n" + pp.stderr.read())
    return pp.returncode

########################################################################################
########################################################################################
########################################################################################

def juage():
    db = pymysql.connect(
        #host="172.16.0.18", 
        host="localhost",
        user="vmtest", 
        password="vmtest888",
        database="weiyuan", 
        autocommit=True,
        charset="utf8")

    #使用 cursor() 方法创建一个游标对象 cursor
    cursor = db.cursor()

    # 分析中状态
    analyse = "ANALYSEING"
    sql_analyse = "select * from tb_analyse_task where status='%s'" %analyse

    cursor.execute(sql_analyse)
    arr_analyse = cursor.fetchone()

    # 排队状态
    queue = "NOT_START"
    sql_quese = "select * from tb_analyse_task where status='%s'" %queue

    cursor.execute(sql_quese)
    arr = cursor.fetchone()


    # 判断是否有任务在运行中，以及是否有任务在排队
    if not arr_analyse: 
        print ("没有任务处于运行中")
        if not arr:
            print ("没有任务在排队")
        else:
            task_id = arr[0]  # 任务id，task_id
            
            # 样本信息表,提取task对应的样本id，以及样本的信息
            sql_task_sample = "select * from tb_task_sample where task_id='%s'" %task_id
            cursor.execute(sql_task_sample)
            cds = cursor.fetchall()
            df_all_sample = pd.DataFrame(list(cds))

            # 第7列为样本id，存在list1上面，后续根据列表来返回html等。
            sample_lists = df_all_sample[7].tolist()
            print (sample_lists)

            mkdir_shell = "mkdir -p /sdbb/Earth/Rawdata/%s /sdbb/Earth/Analysis/%s /sdbb/Earth/AnalysisResults/%s" %(task_id,task_id,task_id) 
            os.system(mkdir_shell)

            # 第7,23,24列分别为name、fq1、fq2，大部分流程只需要这三个数据        
            sample = "/sdbb/Earth/Rawdata/%s/task_sample.csv" %(task_id)
            df_all_sample = df_all_sample[[7,23,24]]
            df_all_sample.to_csv(sample,header=False,index=False,sep='\t')

            # 存放log文件
            log = "/sdbb/Earth/AnalysisResults/%s/log.txt" %task_id
            f = open(log,'w')

            # 运行时，先把运行状态改成分析中
            sql_analyse = "update tb_analyse_task set status='ANALYSEING' where id='%s'" %task_id
            cursor.execute(sql_analyse)
            db.commit()

            sequencing_platform = arr[7]    #测序平台
            sequencing_length = arr[10]  #测序长度
            created_library_strategy = arr[11]  #建库策略
            pathogene = arr[12] # 病原体
            type = arr[13]  #类型
            
            # 鉴定
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            if (type == "APPRAISAL_ANALYSIS"):
                check_cmd = "perl /sdbb/share/pipeline/file/file_check.pl %s > %s" %(sample,log)
                print (check_cmd)
                os.system(check_cmd)
                if (not os.path.getsize(log)):
                    ###############################################################################################################
                    ###############################################################################################################
                    # 鉴定-Nanopore鉴定
                    if (sequencing_platform == "Nanopore"):
                        print ("任务",task_id,time.asctime(),"开始nanopore鉴定流程")
                        shell = "perl /sdbb/share/pipeline/nanopore/bin/nano_identify.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" %(task_id,sample,task_id)
                        print (shell)

                        workflow_cn = "nanopore鉴定"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id)

                    sql_end = "update tb_analyse_task set status='COMPLETE' where id='%s'" %task_id
                    cursor.execute(sql_end)
                    db.commit()
                    
                else:
                    sql_error = "update tb_analyse_task set status='ERROR' where id='%s'" %task_id
                    print (sql_error)
                    cursor.execute(sql_error)
                    db.commit()
                    
                    message = ""
                    with open(log, 'r', encoding='utf-8') as f:  # 打开文件
                        lines = f.readlines()  # 读取所有行
                        #first_line = lines[0].strip()  # 取第一行
                        for line in lines:
                            message += line

                    sql_error = "update tb_analyse_task set error_message='%s' where id='%s'" %(message,task_id)
                    print (sql_error)
                    cursor.execute(sql_error)
                    db.commit()

                    return None
            # 全基因组
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            elif (type == "NOVEL_CORONAVIRUS"):
                check_cmd = "perl /sdbb/share/pipeline/file/file_check.pl %s > %s" %(sample,log)
                print (check_cmd)
                os.system(check_cmd)
                if (not os.path.getsize(log)):
                    ###############################################################################################################
                    ###############################################################################################################
                    # 全基因组-Nanopore新冠
                    if (sequencing_platform == "Nanopore" and pathogene == "2019新型冠状病毒"):
                        print ("任务",task_id,time.asctime(),"开始Nano新冠全基因组流程")
                        shell = "perl /sdbb/share/pipeline/nanopore/bin/nano_sars.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" %(task_id,sample,task_id)
                        print (shell)

                        workflow_cn = "Nano新冠全基因组"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id)
                    ###############################################################################################################
                    ###############################################################################################################
                    # 全基因组-新冠扩增子
                    elif created_library_strategy == "扩增子" and sequencing_platform in ["Illumina", "MGI"] \
                        and pathogene == "2019新型冠状病毒":
                        print ("任务",task_id,time.asctime(),"开始二代测序新冠扩增子分析流程")
                        # bash /sdbb/share/pipeline/Sars_Cov2_Amplicon/sars_cov2_wrap_for_earthFQ.sh task_sample.csv 233
                        shell = "bash /sdbb/share/pipeline/Sars_Cov2_Amplicon/sars_cov2_wrap_for_earthFQ.sh %s %s" \
                            %(sample,task_id)
                        print (shell)
                        
                        workflow_cn = "二代测序新冠扩增子分析"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id, \
                            fa_pattern="/sdbb/Earth/AnalysisResults/%s/%s/4.consensus/%s.consensus.fa")
                    ###############################################################################################################
                    ###############################################################################################################
                    # 单菌、单真菌、单病毒
                    # 全基因组-单细菌
                    elif (created_library_strategy == "单细菌"):
                        
                        print ("任务",task_id,time.asctime(),"开始单细菌denovo组装流程")
                        shell =f"bash /sdbb/share/pipeline/single_denovo/parallel_run.sh {sample} /sdbb/Earth/Analysis/{task_id} bacteria /sdbb/Earth/AnalysisResults/{task_id}"
                        print (shell)

                        workflow_cn = "单细菌denovo组装"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id)
                    ###############################################################################################################
                    # 全基因组-单真菌
                    elif (created_library_strategy == "单真菌"):
                        print ("任务",task_id,time.asctime(),"开始单真菌denovo组装流程")
                        shell =f"bash /sdbb/share/pipeline/single_denovo/parallel_run.sh {sample} /sdbb/Earth/Analysis/{task_id} fungi /sdbb/Earth/AnalysisResults/{task_id}"
                        print (shell)

                        workflow_cn = "单真菌denovo组装"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id)
                    ###############################################################################################################
                    # 全基因组-单病毒
                    elif (created_library_strategy == "单病毒"):
                        print ("任务",task_id,time.asctime(),"开始单病毒denovo组装流程")
                        shell =f"bash /sdbb/share/pipeline/single_denovo/parallel_run.sh {sample} /sdbb/Earth/Analysis/{task_id} virus /sdbb/Earth/AnalysisResults/{task_id}"
                        print (shell)

                        workflow_cn = "单病毒denovo组装"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id)
                    ###############################################################################################################
                    ###############################################################################################################
                    # 全基因组-新冠宏基因组比对-获取全基因组PE150
                    elif (created_library_strategy == "宏基因组" and pathogene == "2019新型冠状病毒" and sequencing_length == "PE150"):
                        print ("任务",task_id,time.asctime(),"开始新冠宏基因组PE150拼接流程") 
                        shell =f"bash /sdbb/share/pipeline/nCov2019_meta_assemble/meta_ncov_2019_pe_assemble/parallel_run.sh {sample} /sdbb/Earth/Analysis/{task_id} /sdbb/Earth/AnalysisResults/{task_id}"
                        print (shell)

                        workflow_cn = "新冠宏基因组PE150拼接"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id)
                    ###############################################################################################################        
                    # 全基因组-新冠宏基因组比对-获取全基因组PE150
                    elif (created_library_strategy == "宏基因组" and pathogene =="2019新型冠状病毒" and sequencing_length == "SE50"):
                        print ("任务",task_id,time.asctime(),"开始新冠宏基因组SE50拼接流程") 
                        shell =f"perl /sdbb/share/pipeline/ncov2019_blast_denovo/bin/ncov2019_blast_denovo.pl -workdir /sdbb/Earth/Analysis/%s -sequence SE -fq_input %s -task_id %s" %(task_id,sample,task_id)
                        print (shell)

                        workflow_cn = "新冠宏基因组SE50拼接"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id) 
                    ###############################################################################################################
                    ###############################################################################################################
                    # 全基因组-宏基因组组装
                    elif (created_library_strategy == "宏基因组" and pathogene =="无参病原"):
                        print ("任务",task_id,time.asctime(),"开始宏基因组组装流程")
                        shell = "perl /sdbb/share/pipeline/meta_denovo/bin/meta_denovo_and_anno.pl -workdir /sdbb/Earth/Analysis/%s -fq_input %s -task_id %s" %(task_id,sample,task_id)
                        print (shell)

                        workflow_cn = "宏基因组组装"
                        popen_returncode = wrap_popen(shell, task_id, f, cursor)
                        if popen_returncode == 0:
                            return2mysql(sample_lists=sample_lists,cursor=cursor,task_id=task_id)     
                
                else:
                    sql_error = "update tb_analyse_task set status='ERROR' where id='%s'" %task_id
                    print (sql_error)
                    cursor.execute(sql_error)
                    db.commit()
                    
                    message = ""
                    with open(log, 'r', encoding='utf-8') as f:  # 打开文件
                        lines = f.readlines()  # 读取所有行
                        #first_line = lines[0].strip()  # 取第一行
                        for line in lines:
                            message += line

                    sql_error = "update tb_analyse_task set error_message='%s' where id='%s'" %(message,task_id)
                    print (sql_error)
                    cursor.execute(sql_error)
                    db.commit()

                    return None

            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            # 16S
            elif (type == "SIXTEEM_S"): 
                check_cmd = "perl /sdbb/share/pipeline/file/file_check.pl %s > %s" %(sample,log)
                print (check_cmd)
                os.system(check_cmd)
                if (not os.path.getsize(log)):

                    print("任务", task_id, time.asctime(), "开始16S测序流程")
                    # Add by renchaobo
                    d_out = f"/sdbb/Earth/Analysis/{task_id}"
                    f_config = config4amplicon(arr, cds, d_out)
                    shell = f"python /sdbb/share/pipeline/16S/main.py -c {f_config} -o {d_out} --run"
                    print (shell)
                    
                    workflow_cn = "16S测序"
                    popen_returncode = wrap_popen(shell, task_id, f, cursor)
                    if popen_returncode == 0:
                        return2mysql_16s(sample_lists=sample_lists,cursor=cursor,task_id=task_id)

                    cmd = f"cp -rf {d_out}/{task_id} /sdbb/Earth/AnalysisResults/{task_id}\ncp {d_out}/{task_id}.zip /sdbb/Earth/AnalysisResults/{task_id}/{task_id}.zip"
                    os.system(cmd)

                else:
                    sql_error = "update tb_analyse_task set status='ERROR' where id='%s'" %task_id
                    print (sql_error)
                    cursor.execute(sql_error)
                    db.commit()
                    
                    message = ""
                    with open(log, 'r', encoding='utf-8') as f:  # 打开文件
                        lines = f.readlines()  # 读取所有行
                        for line in lines:
                            message += line

                    sql_error = "update tb_analyse_task set error_message='%s' where id='%s'" %(message,task_id)
                    print (sql_error)
                    cursor.execute(sql_error)
                    db.commit()

                    return None

            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            # 进化树
            elif (type == "EVOLUTIONARY_TREE"):
                print("任务", task_id, time.asctime(), "开始进化树流程")
                d_analysis = f"/sdbb/Earth/Analysis/{task_id}"
                d_result = f"/sdbb/Earth/AnalysisResults/{task_id}"
                f_fasta = cds[0][36]
                stype = "nucl" if arr[24] == "核酸" else "prot"
                method = arr[22].strip().split('：')[0]
                shell = f"python /sdbb/share/Develop/PhyTreePlugin/main.py -fa {f_fasta} -s {stype} -m {method} -o {d_analysis} >{d_analysis}/run.log 2>&1"

                workflow_cn = "进化树"
                popen_returncode = wrap_popen(shell, task_id, f, cursor)
                os.makedirs(f"{d_result}/Upload", exist_ok=True)
                cmd = f"cp -rfL {d_analysis}/Upload/* {d_result}/Upload\ncp {d_analysis}/result.zip {d_result}/{task_id}.zip"
                os.system(cmd)
                if popen_returncode == 0:
                    return2mysql_16s(sample_lists=sample_lists, cursor=cursor, task_id=task_id,
                        html_pattern=r"/file/weiyuan/AnalysisResults/%s/Upload/result.svg")
            
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            # blast
            elif (type == "BLAST"):
                print("任务", task_id, time.asctime(), "开始Blast流程")
                d_analysis = f"/sdbb/Earth/Analysis/{task_id}"
                d_result = f"/sdbb/Earth/AnalysisResults/{task_id}"
                if os.path.exists(cds[0][37]):
                    f_fasta = cds[0][37]
                else:
                    f_fasta = os.path.join(d_analysis, "sequence.fasta")
                    with open(f_fasta, 'w') as OUT:
                        print(cds[0][37], file=OUT)
                qtype = "nucl" if arr[24] == "核酸" else "prot"
                database = "nr" if arr[17] == "NR" else "nt"
                species = arr[21].strip()
                shell = f"python /sdbb/share/Develop/BlastPlugin/main.py -fa {f_fasta} -t {qtype} -d {database} --species {species} -o {d_analysis}"

                workflow_cn = "Blast"
                popen_returncode = wrap_popen(shell, task_id, f, cursor)
                if popen_returncode == 0:
                    return2mysql_16s(sample_lists=sample_lists, cursor=cursor, task_id=task_id,
                        html_pattern=r"/file/weiyuan/AnalysisResults/%s/Upload/blast.result.tsv")

                os.makedirs(f"{d_result}/Upload", exist_ok=True)
                cmd = f"cp -rfL {d_analysis}/Upload/blast.result.tsv {d_result}/Upload;cp -rf {d_analysis}/result.zip {d_result}/{task_id}.zip"
                os.system(cmd)

            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            else:
                print ("任务",task_id,time.asctime(),"未知流程")
            
            # 最后
            if popen_returncode is None:
                # 终止
                f.write("任务%s %s 终止%s流程" %(task_id,time.asctime(), workflow_cn))
                sql_end = "update tb_analyse_task set status='CANCEL' where id='%s'" %task_id
            elif popen_returncode == 0:
                # 完成
                sql_end = "update tb_analyse_task set status='COMPLETE' where id='%s'" %task_id
                f.write("任务%s %s 完成%s流程" %(task_id,time.asctime(), workflow_cn))
            else:
                # 异常中断
                sql_end = "update tb_analyse_task set status='ERROR' where id='%s'" %task_id
                f.write("任务%s %s 异常中断%s流程" %(task_id,time.asctime(), workflow_cn))

            cursor.execute(sql_end)
            db.commit()
            

    else:
        print ("有任务在运行中")

    cursor.close()  # 关闭游标
    db.close()  # 关闭数据库连接


def loop_func(func, second):
 # 每隔second秒执行func函数
    while True:
        func()
        time.sleep(second)


loop_func(juage,10)

