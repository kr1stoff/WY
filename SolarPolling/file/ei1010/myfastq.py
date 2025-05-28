import sys
import logging
import gzip
from pathlib import Path
# custom
from file.ei1010 import common


# 给库添加日志
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class AnyTo100BP():
    def __init__(self, origin_fq, trimmed_fq):
        """
        [20230221 mxf] 任意长度FASTQ修剪100bp, 使用'seqtk trimfq'
        1. 如果小于100则保留
        2. 如果长度不统一按照最长的减去100的差值来修剪
        3. 无论输入FASTQ是否压缩, 输出都不压缩.fastq后缀
        """
        self.origin_fq = origin_fq
        self.trimmed_fq = trimmed_fq

    def guess_suffix_open(self):
        #猜测压缩/未压缩的FASTQ
        self.guess_open = gzip.open if Path(self.origin_fq).suffix == '.gz' else open

    def guess_read_length(self):
        #猜测读长, fastq文件前40行猜测测序长度
        line_number = 0
        lengths = []
        with self.guess_open(self.origin_fq, 'rt') as f:
            for line in f.readlines():
                line_number += 1
                if (line_number % 4) in [1,3]: continue # 1,3行不是序列长度
                if (line_number % 4) in [0,2]: lengths.append(len(line.strip()))
                if line_number >= 40: break
        self.seq_len = max(lengths)

    def any_to_100bp_command(self):
        #执行修剪
        params = common.get_params()
        self.guess_suffix_open()
        self.guess_read_length()
        trim_length = self.seq_len - 101 # 按照最长的read长度减101bp算修剪长度
        if trim_length > 1:
            if self.trimmed_fq.endswith('gz'):
                cml = f'{params.seqtk} trimfq -e {trim_length} {self.origin_fq} | gzip > {self.trimmed_fq}'
            else:
                cml = f'{params.seqtk} trimfq -e {trim_length} {self.origin_fq} > {self.trimmed_fq}'
        else:
            cml = f'ln -sf {self.origin_fq} {self.trimmed_fq}'
        return cml
