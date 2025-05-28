# README

微远公卫分型溯源流程

*添加物种需要修改文件`config/content.tsv`,数据库表`tb_dict_genome_type`、`tb_dict_ref_gene`*

## :peach: 环境依赖

- R4.1.2.tar.gz
- snippy.tar.gz
- mlva.tar.gz
- hicap.tar.gz
- typetrace.tar.gz
- shigatpyer.tar.gz

## :rice: 分型

### MLST

MLST分型

### SeroType

血清型分析

### Trace

### MLVA

序列文件应以fa结尾

```bash
python mlva.py -i {d_in} -o {d_out} -p {f_primer}
```

### 其它分型方法

#### 金黄色葡萄球菌的Spa分型

*可直接使用脚本，也可按照模块依照需求导入*
[脚本逻辑来源](https://github.com/zwets/spa-type)

```bash
python lib/spa.py -n test -f t.fa --primer /sdbb/bioinfor/renchaobo/test/t/spa-type/reference/primers.fna --spacsv /sdbb/bioinfor/renchaobo/test/t/spa-type/database/spa_db.csv --spatype /sdbb/bioinfor/renchaobo/test/t/spa-type/database/spatypes.txt -o spa
```

溯源分析

## :bread: 其它功能

### Enterobase数据库下载与构建

列出所有可用的scheme

```bash
python EnteroBase.py list
```

下载目标scheme

```bash
python EnteroBase.py download --scheme Escherichia.cgMLSTv1 -o Escherichia.cgMLSTv1
```

scheme参数为Enterobase可用的scheme

构建cgMLST数据库

```bash
python EnteroBase.py build -s profiles.list.gz -i Streptococcus$i -o Streptococcus -p Streptococcus
```

### Pasteur数据库下载与构建

列出所有可用的scheme

```bash
python Pasteur.py list
```

将会下载Pasteur数据库中API可用的scheme的信息到scheme.json文件中，由于改数据库API的问题，会导致一些在网页端可见的数据库在scheme.json文件中并不可见，
需要用户切换到物种具体的网页获取数据库名称，类似 https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_diphtheria_seqdef中的数据库名称pubmlst_diphtheria_seqdef，
通过构建API地址 https://bigsdb.pasteur.fr/api/db/pubmlst_diphtheria_seqdef/schemes 获取可用的schemes及相关信息

```bash
python Pasteur.py download 
```

### :cry: xSNP数据库下载与构建

1. 下载目标物种的基因组序列和注释文件

```bash
# target.list为目标物种拉丁名，一行为一个物种
cat target.list | while read TARGET; do
  name=$(echo "$TARGET" | tr ' ' '_')
  mkdir -p $name
  awk -F'\t' -v TARGET="$TARGET" '($5=="reference genome" || $5=="representative genome") && $8 ~ TARGET' /sdbb/bioinfor/renchaobo/Project/WY2022111601/Work/00.Prepare/assembly_summary_refseq.txt >$name/info.tsv
done
for i in */info.tsv; do (
  cd $(dirname $i)
  cat info.tsv | while read line; do
    path=$(echo "$line" | awk -F'\t' '{print $20}')
    name=$(basename $path)
    dlink="$(echo $path | sed s,https://ftp.ncbi.nlm.nih.gov,anonftp@ftp.ncbi.nlm.nih.gov:,g)/${name}_genomic.fna.gz"
    /home/renchaobo/.aspera/connect/bin/ascp -i /home/renchaobo/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l100m $dlink ./
  done
); done
for i in */info.tsv; do (
  cd $(dirname $i)
  cat info.tsv | while read line; do
    path=$(echo "$line" | awk -F'\t' '{print $20}')
    name=$(basename $path)
    dlink="$(echo $path | sed s,https://ftp.ncbi.nlm.nih.gov,anonftp@ftp.ncbi.nlm.nih.gov:,g)/${name}_genomic.gff.gz"
    /home/renchaobo/.aspera/connect/bin/ascp -i /home/renchaobo/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l100m $dlink ./
  done
); done
```

2. 将基因组序列解压并重命名为ref.fna，将注释文件解压并提取gene信息，转换为bed格式，生成mask.bed文件

```bash
for i in /sdbb/bioinfor/renchaobo/Project/WY2024010501/Work/02.Database/*/*.fna.gz;do (dname=$(dirname $i);sname=$(basename $dname);mkdir -p $sname;gunzip -c $i > /sdbb/bioinfor/renchaobo/Project/WY2024010501/Work/03.xSNP/$sname/ref.fna);done
for i in /sdbb/bioinfor/renchaobo/Project/WY2024010501/Work/02.Database/*/*.gff.gz;do (dname=$(dirname $i);sname=$(basename $dname);zcat $i|awk -F'\t' -vOFS="\t" '$3=="gene"' -|gff2bed - > /sdbb/bioinfor/renchaobo/Project/WY2024010501/Work/03.xSNP/$sname/target.bed);done
for i in */*.fna;do python /sdbb/bioinfor/renchaobo/Develop/ToolBox/fa/SeqTools.py length -f $i -o ${i%.*}.info;done
# 部分物种需要修改ref.info文件内染色体顺序
for i in */*.info;do dname=$(dirname $i);(cd $dname;bedtools complement -g ref.info -i target.bed > mask.bed);done
```

