# 二代测序扩增子通用分析流程
## 用法
```
python AmpliconGP/main.py -i AmpliconGP/template/input.yaml
```

## 更新
- [231229] 流程更新
  - microbial_wgs 环境要更新 (click)
  - fastqc, bam, consensus_fa 展示
  - dos2unix ref.fa, 换行符问题优化 
- [221031] 流程优化  
    过滤前后fastqc, 添加段落icon  
- [221012] 流程逻辑优化  
    1)类似新冠,把fastp参数拿出来; 2)FQ质控表格统一; 3)流程图优化
