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

    def __init__(self, ptype: str):
        """
        Parser for Primer3 result
        """
        self.ptype = ptype
        self.header = self.get_header(self.ptype)
        logger.debug(self.header)

    def get_header(self, ptype: str):
        """
        获取特定类型引物设计输出表格的表头

        :param ptype: 引物设计类型[ internal | external]
        """
        if ptype == "internal":
            header = ["Index", 
                      "ForwardStart", "ForwardEnd",
                      "ForwardLength",
                      "ForwardTM", "ForwardGC", "ForwardSeq",
                      "ReverseStart", "ReverseEnd",
                      "ReverseLength",
                      "ReverseTM", "ReverseGC", "ReverseSeq",
                      "ProbeStart", "ProbeEnd", "ProbeLength",
                      "ProbeTM", "ProbeGC", "ProbeSeq",
                      "AmpliconTM", "AmpliconGC", "AmpliconTa", "AmpliconLength", "PenaltyPrimerExpress", 
                      "Chromosome", "Start", "End"]
        elif ptype == "external":
            header = ["Index", "Chromosome",
                      "ForwardStart", "ForwardEnd",
                      "ForwardLength",
                      "ForwardTM", "ForwardGC", "ForwardSeq",
                      "ReverseStart", "ReverseEnd",
                      "ReverseLength",
                      "ReverseTM", "ReverseGC", "ReverseSeq"]
        else:
            sys.exit(logger.error(f"Unsupport type {ptype} for primer3 result"))
        return header

    def parse_sequence_id(self, seqid):
        """解析sequence_id, 符合格式`chrom:start-end`则返回[chrom,start,end]"""
        ptn = re.compile(r'(.*?):(\d+)-(\d+)')
        if re.match(ptn, seqid):
            return re.findall(ptn, seqid)[0]
        else:
            return (seqid, '-', '-')

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
        amplicon_tm = re.findall(r"PRIMER_PAIR_\d+_PRODUCT_TM=(\S+)", content)
        amplicon_gc = '-'
        amplicon_ta = '-'
        if self.ptype == "internal":
            probe_pos = re.findall(r"PRIMER_INTERNAL_\d+=(\d+),(\d+)", content)
            probe_seq = re.findall(r"PRIMER_INTERNAL_\d+_SEQUENCE=(\S+)", content)
            probe_tm = re.findall(r"PRIMER_INTERNAL_\d+_TM=(\S+)", content)
            probe_gc = re.findall(r"PRIMER_INTERNAL_\d+_GC_PERCENT=(\S+)", content)

        tmp = []
        for index in range(len(forward_pos)):
            if self.ptype == "internal":
                #PrimerExpress罚分. 最优引物长度20,罚分系数1.0; 最小扩增子长度70,罚分系数5.
                opt_prm_len, len_mult, min_amp_len, min_amp_mult = 20, 1.0, 70, 5
                penalty_prmexp = (abs((int(forward_pos[index][1]) - opt_prm_len) * len_mult)) + \
                (abs((int(reverse_pos[index][1]) - opt_prm_len) * len_mult)) + \
                ((int(amplicon_len[index]) - min_amp_len) * min_amp_mult)
                #解析sequence_id
                chrom, start, end = self.parse_sequence_id(sequence_id)
                res = [index, 
                       int(forward_pos[index][0]),
                       int(forward_pos[index][0]) + int(forward_pos[index][1]),
                       int(forward_pos[index][1]), 
                       forward_tm[index], forward_gc[index], forward_seq[index],
                       int(reverse_pos[index][0]) + 1, #Reverse Start大 End小
                       int(reverse_pos[index][0]) - int(reverse_pos[index][1]) + 1,
                       int(reverse_pos[index][1]), 
                       reverse_tm[index], reverse_gc[index], reverse_seq[index],
                       int(probe_pos[index][0]),
                       int(probe_pos[index][0]) + int(probe_pos[index][1]),
                       int(probe_pos[index][1]), 
                       probe_tm[index], probe_gc[index], probe_seq[index],
                       amplicon_tm[index], amplicon_gc, amplicon_ta, amplicon_len[index], penalty_prmexp,
                       chrom, start, end]
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
                    print(f">{i}_{info[0]}-forward\n{info[6]}", file=OUT)
                    print(f">{i}_{info[0]}-reverse\n{info[12]}", file=OUT)
                    if self.ptype == "internal":
                        print(f">{i}_{info[0]}-probe\n{info[18]}", file=OUT)

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
                    print(*[i, info[1], info[2], f"{i}_{info[0]}-forward",
                            0, "+"], sep="\t", file=OUT)
                    print(*[i, info[7], info[8], f"{i}_{info[0]}-reverse",
                            0, "-"], sep="\t", file=OUT)
                    if self.ptype == "internal":
                        print(*[i, info[13], info[14],
                                f"{i}_{info[0]}-probe", 0, "+"], sep="\t",
                              file=OUT)
