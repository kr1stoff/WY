# TracePathogen

## 介绍
溯源进化树分析流程，进化溯源结果报告
流程分为三条路线：  
- 核心基因 
- 全基因组SNP
- 全基因组多序列比对  

## 使用方法
### 获取帮助
```
python main.py -h
```
全局参数:  
```
-p {wgs,snp,core}, --pipe {wgs,snp,core}    选择要跑的子流程
-i INYAML, --inyaml INYAML                  样本信息YAML文件.
-n, --dryrun                                不运行,只生成脚本和目录. (default:False)
```

### 应用示例
```
python main.py -p wgs -i template/core.yaml
python main.py -p snp -i template/snpfa.yaml
python main.py -p snp -i template/snpfq.yaml
python main.py -p core -i template/core.yaml
```

## 更新
- [231031] 参数调整 
    `vcf2phylip.py -m MIN_SAMPLES_LOCUS` 调整为样本数 * 0.1
- [x] 重写整个程序，把各个子流程拆分出来
- [x] 全基因组进化树，一个文件一个样本（之前修改过）
- [x] 输出 .tre 或 .nwk (本来就有)
- [221021] 流程升级  
    重写完成
- [220708] 流程优化   
    1.优化细菌SNP流程.加入 `--fasta` 参数, `fasttree` 用 .`fasta` 建树; 2.书写规范优化
- [220616] snp  
    细菌真菌样本SNP特别多, `vcf2phylip.py` 的 `-m MIN_SAMPLES_LOCUS` 参数需要调整,设置为样本数量大小
- [220517] snp  
    流程支持 `FASTQ`，使用与前三种不同的逻辑  
