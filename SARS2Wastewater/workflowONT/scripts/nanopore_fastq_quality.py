#!/usr/bin/env python


import gzip
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def wrap_open(fq_file:str):
    if fq_file.endswith(".gz"):
        whandle = gzip.open(fq_file, "rt")
    else:
        whandle = open(fq_file, "rt")
    return whandle

# fastq[.gz] 文字统计和生成表格
def fq2stats(fq_file, out_read_stats_file, phred=33):
    outlists = list()
    fh = wrap_open(fq_file)
    with open(out_read_stats_file, "w", newline="", encoding="utf-8") as gh:
        header = ["read_id", "read_length", "read_quality"]
        gh.write("\t".join(header) + "\n")
        # 行数记录，每四行是否完整检查
        line_number = 0
        complete_check_list = list()
        for line in fh:
            line_number += 1
            # 逐行看
            if line_number % 4 == 1 and line.startswith("@"):
                complete_check_list.append(0)
                read_id = line.strip().split(" ")[0].replace("@", "")
            elif line_number % 4 == 2:
                complete_check_list.append(1)
            elif line_number % 4 == 3 and line.strip() == "+":
                complete_check_list.append(2)
            elif line_number % 4 == 0:
                complete_check_list.append(3)
                qual_list           = list(line.strip())
                # nanopore Q值 计算方法
                quals = np.fromiter(
                    (ord(x) - int(phred) for x in qual_list),
                    dtype=int, count=len(qual_list))
                mean_read_quality = -10*np.log10(np.mean(10**(quals/-10)))
                read_len = len(qual_list)
            # 每四行是否完整检查
            if len(complete_check_list) == 4:
                outlist = [read_id, read_len, mean_read_quality]
                outlists.append(outlist)
                complete_check_list = list()
                gh.write("\t".join(map(str, outlist)) + "\n")
    fh.close()
    # read 长度,质量表格
    stats_df = pd.DataFrame(outlists, columns=header)
    stats_df.read_quality = stats_df.read_quality.astype(int)
    return stats_df

# 密度图 
def qc_chart(stats_df, out_read_stats_png):
    # seaborn
    sns.set_style("darkgrid")
    ax = sns.jointplot(data=stats_df, x="read_length", y="read_quality",
                kind="kde", fill=True, color="#87CEEB", space=0, ratio=10,
                joint_kws=dict(alpha =0.8,fill=True),
                marginal_kws=dict(fill=True))
    ax.set_axis_labels("Read Length", "Mean Read Quality")
    plt.savefig(out_read_stats_png, dpi=500)
    plt.close()


if __name__ == "__main__":
    fq_file = snakemake.input[0]
    out_read_stats_file = snakemake.output[0]
    out_read_stats_png = snakemake.output[1]
    stats_df = fq2stats(fq_file, out_read_stats_file)
    qc_chart(stats_df, out_read_stats_png)
