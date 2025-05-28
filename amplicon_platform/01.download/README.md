# README

微远基因溯源平台

## 预扩增引物设计
内外引物及探针设计

----
使用配置文件篇`pPCR.yml`，现在只有内引物及探针的设计参数有用，外引物设计采用PrimerPlex
运行命令

```bash
python main.py ppcr -c pPCR.yml -o ./
```

----
TODO ：

1. 外引物设计所使用的PrimerPlex,如果无法在目标区域周围找到引物，则会报错，需要修改PrimerPlex或者使用Primer3重写
2. SNP突变检测现阶段显得可有可无

## 探针设计

None