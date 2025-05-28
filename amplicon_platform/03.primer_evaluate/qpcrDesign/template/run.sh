cd /sdbb/bioinfor/mengxf/TASKS/WY23051801
#单个
python /sdbb/bioinfor/mengxf/Software/gitee/WY/qpcrDesign/main.py qpcr -f /sdbb/bioinfor/mengxf/TASKS/WY23051801/446/region.part_010/Prepare/sequence.fasta -o test
#批量
python /sdbb/bioinfor/mengxf/Software/gitee/WY/qpcrDesign/main.py qpcr-batch -r /sdbb/bioinfor/puzl/Project/WY2022120601/Work/02.Filter/446/ref.fna -b bed/446.result.filter.sort.bed -o 446
