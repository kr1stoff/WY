from encodings import utf_8_sig
import pandas as pd
import os
import sys

IN=os.path.abspath(sys.argv[1])
OUT=os.path.abspath(sys.argv[2])


data=pd.DataFrame(pd.read_csv(IN,encoding='utf_8_sig',sep="\t"))

df=data[["REGION","POS","REF","ALT","TOTAL_DP"]]
df.columns=["参考序列","变异位点","参考序列原碱基","变异后的碱基","突变位点测序深度"]

df.to_excel(f"{OUT}/variants.xlsx",index=False,encoding='utf_8_sig')
df.to_csv(f"{OUT}/variants.txt",sep="\t",index=False,encoding='utf_8_sig')
