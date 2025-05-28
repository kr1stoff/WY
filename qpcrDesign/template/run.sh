# 1. qPCR 引物设计
cd /sdbb/bioinfor/mengxf/TASKS/WY23051801
# 1.1. 单个
python /sdbb/bioinfor/mengxf/Software/gitee/WY/qpcrDesign/main.py qpcr \
    -f /sdbb/bioinfor/mengxf/TASKS/WY23051801/446/region.part_010/Prepare/sequence.fasta \
    -o test
# 1.2. 批量
python /sdbb/bioinfor/mengxf/Software/gitee/WY/qpcrDesign/main.py qpcr-batch \
    -r /sdbb/bioinfor/puzl/Project/WY2022120601/Work/02.Filter/446/ref.fna \
    -b bed/446.result.filter.sort.bed \
    -o 446

# 2. 引物简并分析
cd /sdbb/bioinfor/mengxf/TASKS/WY23091501/440266
python /sdbb/bioinfor/mengxf/Software/gitee/WY/qpcrDesign/main.py translate-degenerate-base \
    -p wu多瘤病毒_merge_primer.xls \
    -a /mnt/VolB02/lanlei/01.Amplicon/01.genome_database_merge/440266/all.fna
