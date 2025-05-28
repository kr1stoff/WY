#!/usr/bin/env python

import os
import sys
import logging
import vcf


# 设置运行日志
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(asctime)s - %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S"
)

# QD < 2.0 || FS > 100.0 || MQ < 40.0 || ReadPosRankSum < -8.0  # SOR > 4.0
def my_snp_filter(record):
    if ("QD" in record.INFO and record.INFO["QD"] < 2.0) \
        or ("FS" in record.INFO and record.INFO["FS"] > 100.0) \
        or ("MQ" in record.INFO and record.INFO["MQ"] < 40.0) \
        or ("ReadPosRankSum" in record.INFO and record.INFO["ReadPosRankSum"] < -8.0):
        record.FILTER = "my_snp_filter"
    else:
        record.FILTER = "PASS"

# DP < 20.0 || QD < 2.0 || FS > 200.0 || SOR > 10.0
def my_indel_filter(record):
    if ("DP" in record.INFO and record.INFO["DP"] < 20)  \
        or ("QD" in record.INFO and record.INFO["QD"] < 2.0) \
        or ("FS" in record.INFO and record.INFO["FS"] > 200.0) \
        or ("SOR" in record.INFO and record.INFO["SOR"] > 10.0):
        record.FILTER = "my_indel_filter"
    else:
        record.FILTER = "PASS"

def main(invcf, filtered_vcf, final_vcf):
    vcf_reader = vcf.Reader(open(invcf, "rt"))
    vcf_writer_filtered = vcf.Writer(open(filtered_vcf,"wt",newline="",encoding="utf-8"), vcf_reader)
    vcf_writer_final = vcf.Writer(open(final_vcf,"wt",newline="",encoding="utf-8"), vcf_reader)
    for record in vcf_reader:
        if record.is_snp: # snp
            my_snp_filter(record)
            vcf_writer_filtered.write_record(record)
        elif record.is_indel: # indel
            my_indel_filter(record)
            vcf_writer_filtered.write_record(record)
        else:
            logging.info(f"其他变异类型: {record} {record.INFO}")
        if record.FILTER == "PASS": # 仅保留PASS的变异
            vcf_writer_final.write_record(record)
    vcf_writer_filtered.close()
    vcf_writer_final.close()

if __name__ == "__main__":
    # 获取参数
    # 输入: gatk原始vcf
    # 输出: 1标记FILTER后的vcf, 2仅保留PASS的vcf
    usage = f"Usage:\n\tpython {os.path.basename(__file__)} <gatk_vcf> <filtered_vcf> <final_vcf>\n"
    if (len(sys.argv) != 4) or (not os.path.isfile(sys.argv[1])):
        print(usage)
        sys.exit(1)
    else:
        invcf, filtered_vcf, final_vcf = sys.argv[1:]
        main(invcf, filtered_vcf, final_vcf)
