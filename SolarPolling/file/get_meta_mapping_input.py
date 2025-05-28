import subprocess
import yaml
import os

def config(arr, cds, d_analysis):
    df[['sample_number','sequencing_data1','sequencing_data2','label','label_compare']]
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
