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

## 一些问题
1. `-region`BED区域在模板靠两边可能会引起`networkx`包的报错, 需要手动剔除部分BED记录.

## 更新
