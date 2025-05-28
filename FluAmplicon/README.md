# FluAmplicon

VENUS 流感扩增子流程

# 可能的问题

1. PE没有组装出来的序列
   IRMA的组装序列的过程中默认会将双端的reads merge到一起，须在配置文件中添加`NO_MERGE=1`
2. SE没有组装出来的序列
   可能需要调整IRMA配置文件的`MIN_LEN`参数


