import subprocess
import yaml
import os

def linkdata(info_sample, d_rawdata):
    """
    Link the rawdata to the analysis dir and rename it to the standard format

    :param info_sample: The sample info table
    :param d_rawdata: The rawdata dir in the pipeline
    """
    samples = []
    for record in info_sample:
        sample_name = record['sample_number']
        if (record['sequencing_data2'] != ""):
            cmd = f"""ln -sf {record['sequencing_data1']} {d_rawdata}/{sample_name}_1.fq.gz
ln -sf {record['sequencing_data2']} {d_rawdata}/{sample_name}_2.fq.gz"""
            subprocess.run(cmd, shell=True)
        elif (record['sequencing_data2'] == ""):
            cmd = f"ln -sf {record['sequencing_data1']} {d_rawdata}/{sample_name}_1.fq.gz"
            subprocess.run(cmd, shell=True)
        else:
            raise ValueError(f"UNSUPPORT SAMPLETYPE")
        samples.append(sample_name)
    return samples


def get_group_info(info_sample):
    """
    Get the group info

    :param info_sample: The sample info table
    """
    res = {}
    for record in info_sample:
        sample_name = record['sample_number']
        group_name = record['label']
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
        info_compare = record['label_compare']
        if info_compare:
            for compare in info_compare.strip().split(','):
                res[tuple(compare)] = 1
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
    print (info_task['task_id'])
    res["project"] = info_task['task_id']
    # 测序区域
    region_map = {"ITS": "ITS2", "18S": "V4","":""}
    info_type = info_task['rrna_type'].strip().split(',')
    print (info_type)
    res["type"] = info_type[0]
    res["region"] = info_type[1] if info_type[0] == "16S" else region_map[info_type[0]]
    # 数据库
    db_map = {"16S": "greengenes", "18S": "silva", "ITS": "unite", "": ""}
    res["database"] = db_map[info_type[0]]
    # 样本路径
    d_data = os.path.join(d_analysis, "rawdata")
    os.makedirs(d_data, exist_ok=True)
    res["rawdata"] = d_data
    res["samples"] = linkdata(info_sample, d_data)
    # 数据类型
    res["data_type"] = "PE" if info_task['pe_or_se'].startswith("双端") else "SE"
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
