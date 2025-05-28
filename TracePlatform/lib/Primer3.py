#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/2/14 10:28
# @Last Modified by:   Ming
# @Last Modified time: 2022/2/14 10:28
import logging
import re
import sys
from pathlib import Path

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

d_bin = Path(__file__).absolute().parent


class Setting():
    """
    Primer3 setting object
    """

    def __init__(self, name: str, sequence: str):
        """
        Primer3 setting

        :param name: The sequence name
        :param sequence: The fasta sequence to fetch primer
        """
        # 项目配置
        self.task_conf = {}
        self.set("SEQUENCE_ID", name)
        self.set("SEQUENCE_TEMPLATE", sequence)

    def __str__(self):
        """
        For print the content
        """
        res = []
        for key, value in self.task_conf.items():
            res.append(f"{key}={value}")
        res.append("=")
        return "\n".join(res)

    def set(self, key: str, value: str):
        """
        Set the paramter for Primer3

        :param key: The parameter for Primer3
        :param value: The value for Primer3
        """
        self.task_conf[key] = value

    def write(self, f_out: Path):
        """
        Out put the Primer3 setting to file

        :param f_out: The file name of output
        """
        with open(f_out, 'w') as OUT:
            for key, value in self.task_conf.items():
                print(f"{key}={value}", file=OUT)
            print("=", file=OUT)


class Parser():
    """
    Primer3 结果解析器

    :TODO 不是很灵活，后续更改
    """

    def __init__(self, ptype: str, 
                 opt_prm_len:int, len_mult:float, 
                 min_amp_len:int, min_amp_mult:int):
        """
        Parser for Primer3 result
        """
        self.ptype = ptype
        self.header = self.get_header(self.ptype)
        logger.debug(self.header)
        self.opt_prm_len = opt_prm_len
        self.len_mult = len_mult
        self.min_amp_len = min_amp_len
        self.min_amp_mult = min_amp_mult

    def get_header(self, ptype: str):
        """
        获取特定类型引物设计输出表格的表头

        :param ptype: 引物设计类型[ internal | external]
        """
        if ptype == "internal":
            header = ["Index", "Chromosome",
                      "ForwardStart(0base)", "ForwardEnd(1base)",
                      "ForwardLength",
                      "ForwardTM", "ForwardGC", "ForwardSeq",
                      "ReverseStart(0base)", "ReverseEnd(1base)",
                      "ReverseLength",
                      "ReverseTM", "ReverseGC", "ReverseSeq",
                      "ProbeStart(0base)", "ProbeEnd(1base)", "ProbeLength",
                      "ProbeTM", "ProbeGC", "ProbeSeq",
                      "AmpliconLength", "PenaltyPrimerExpress"]
        elif ptype == "external":
            header = ["Index", "Chromosome",
                      "ForwardStart(0base)", "ForwardEnd(1base)",
                      "ForwardLength",
                      "ForwardTM", "ForwardGC", "ForwardSeq",
                      "ReverseStart(0base)", "ReverseEnd(1base)",
                      "ReverseLength",
                      "ReverseTM", "ReverseGC", "ReverseSeq"]
        else:
            sys.exit(logger.error(f"Unsupport type {ptype} for primer3 result"))
        return header

    def load_file(self, f_name: Path):
        """
        Load the Primer3 result file

        :param f_name: The Primer3 result file
        """
        self.info = {}

        content = open(f_name, 'r').read()
        sequence_id = re.findall(r"SEQUENCE_ID=(\S+)", content)[0]
        forward_pos = re.findall(r"PRIMER_LEFT_\d+=(\d+),(\d+)", content)
        forward_seq = re.findall(r"PRIMER_LEFT_\d+_SEQUENCE=(\S+)", content)
        forward_tm = re.findall(r"PRIMER_LEFT_\d+_TM=(\S+)", content)
        forward_gc = re.findall(r"PRIMER_LEFT_\d+_GC_PERCENT=(\S+)", content)
        reverse_pos = re.findall(r"PRIMER_RIGHT_\d+=(\d+),(\d+)", content)
        reverse_seq = re.findall(r"PRIMER_RIGHT_\d+_SEQUENCE=(\S+)", content)
        reverse_tm = re.findall(r"PRIMER_RIGHT_\d+_TM=(\S+)", content)
        reverse_gc = re.findall(r"PRIMER_RIGHT_\d+_GC_PERCENT=(\S+)", content)
        amplicon_len = re.findall(r"PRIMER_PAIR_\d+_PRODUCT_SIZE=(\S+)", content)
        if self.ptype == "internal":
            probe_pos = re.findall(r"PRIMER_INTERNAL_\d+=(\d+),(\d+)", content)
            probe_seq = re.findall(r"PRIMER_INTERNAL_\d+_SEQUENCE=(\S+)",
                                   content)
            probe_tm = re.findall(r"PRIMER_INTERNAL_\d+_TM=(\S+)", content)
            probe_gc = re.findall(r"PRIMER_INTERNAL_\d+_GC_PERCENT=(\S+)",
                                  content)

        tmp = []
        for index in range(len(forward_pos)):
            if self.ptype == "internal":
                # [mengxf 20230511] PrimerExpress罚分
                penalty_prmexp = (abs((int(forward_pos[index][1]) - self.opt_prm_len) * self.len_mult)) + \
                (abs((int(reverse_pos[index][1]) - self.opt_prm_len) * self.len_mult)) + \
                ((int(amplicon_len[index]) - self.min_amp_len) * self.min_amp_mult)
                res = [index, sequence_id,
                       int(forward_pos[index][0]),
                       int(forward_pos[index][0]) + int(forward_pos[index][1]),
                       int(forward_pos[index][1]), forward_tm[index],
                       forward_gc[index], forward_seq[index],
                       int(reverse_pos[index][0]) - int(
                           reverse_pos[index][1]) + 1,
                       int(reverse_pos[index][0]) + 1,
                       int(reverse_pos[index][1]), reverse_tm[index],
                       reverse_gc[index], reverse_seq[index],
                       int(probe_pos[index][0]),
                       int(probe_pos[index][0]) + int(probe_pos[index][1]),
                       int(probe_pos[index][1]), probe_tm[index],
                       probe_gc[index], probe_seq[index],
                       amplicon_len[index], penalty_prmexp]
            elif self.ptype == "external":
                res = [index, sequence_id,
                       int(forward_pos[index][0]),
                       int(forward_pos[index][0]) + int(forward_pos[index][1]),
                       int(forward_pos[index][1]), forward_tm[index],
                       forward_gc[index], forward_seq[index],
                       int(reverse_pos[index][0]) - int(reverse_pos[index][1]),
                       int(reverse_pos[index][0]) + 1,
                       int(reverse_pos[index][1]), reverse_tm[index],
                       reverse_gc[index], reverse_seq[index]]
            else:
                sys.exit(logger.error("Unsupport type for primer3 result"))

            tmp.append(res)
        self.info[sequence_id] = tmp

    def filter(self):
        """
        引物过滤

        :TODO 添加过滤参数
        """

        def is_good_primer(sequence: str):
            """
            判断引物是否合格

            :param sequence: The primer sequence
            """
            flag = 0
            for i in sequence[:-5]:
                if i.islower():
                    flag += 1
            if flag > 2:
                return False

            for i in sequence[-5:]:
                if i.islower():
                    return False
            return True

        def is_good_probe(sequence: str):
            """
            判断探针是否合格

            :param sequence: The probe sequence
            """
            for i in sequence[:5]:
                if i.islower():
                    return False
            for i in sequence[-5:]:
                if i.islower():
                    return False
            snp_num = 0
            for i in sequence[5:-5]:
                if i.islower():
                    snp_num += 1
            if snp_num > 1:
                return False
            return True

        tmp = {}
        for i, j in self.info.items():
            tmp[i] = []
            for info in j:
                if is_good_primer(info[7]) and is_good_primer(info[13]):
                    if self.ptype != "internal":
                        tmp[i].append(info)
                    if self.ptype == "internal" and is_good_probe(info[19]):
                        tmp[i].append(info)
        self.info = tmp

    def filter_region(self, start: int, end: int):
        """
        Remove the region that the start is small than left start and bigger than right end
        """
        tmp = {}
        for i, j in self.info.items():
            tmp[i] = []
            for info in j:
                if info[2] >= start and info[9] <= end:
                    tmp[i].append(info)
        self.info = tmp

    def to_csv(self, f_out: Path, sep="\t"):
        """
        Write all of Primer3 result to a tsv file

        :param f_out: The output file name
        :param seq: The output delimiter[default: '\t']
        """
        logger.info(f"Output the primer3 result info to {f_out}")
        with open(f_out, 'w') as OUT:
            print(*self.header, sep=sep, file=OUT)
            for i, j in self.info.items():
                for info in j:
                    print(*info, sep=sep, file=OUT)

    def to_fasta(self, f_out: Path):
        """
        Write the primer to fasta file

        :param f_out: The output fasta file name
        """
        logger.info(f"Out put the fasta info to {f_out}")
        with open(f_out, 'w') as OUT:
            for i, j in self.info.items():
                for info in j:
                    logger.debug(info)
                    print(f">{i}_{info[0]}-forward\n{info[7]}", file=OUT)
                    print(f">{i}_{info[0]}-reverse\n{info[13]}", file=OUT)
                    if self.ptype == "internal":
                        print(f">{i}_{info[0]}-probe\n{info[19]}", file=OUT)

    def to_bed(self, f_out: Path):
        """
        Write the primer to bed file

        :param f_out: The output bed file name
        """
        logger.info(f"Output the bed info to {f_out}")
        with open(f_out, 'w') as OUT:
            for i, j in self.info.items():
                for info in j:
                    logger.debug(info)
                    print(*[info[1], info[2], info[3], f"{i}_{info[0]}-forward",
                            0, "+"], sep="\t", file=OUT)
                    print(*[info[1], info[8], info[9], f"{i}_{info[0]}-reverse",
                            0, "-"], sep="\t", file=OUT)
                    if self.ptype == "internal":
                        print(*[info[1], info[14], info[15],
                                f"{i}_{info[0]}-probe", 0, "+"], sep="\t",
                              file=OUT)
