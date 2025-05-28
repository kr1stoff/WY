import sys
import os
import pandas as pd

infile=os.path.abspath(sys.argv[1])
depfile=os.path.abspath(sys.argv[2])
outdir=os.path.dirname(infile)


if os.path.exists(infile) and os.path.exists(depfile):

    col=["Contig","ARO","Best_Hit_Bitscore","Best_Identities","Percentage Length of Reference Sequence","Drug Class","Resistance Mechanism",]
    #选取部分列
    rgi_df=pd.DataFrame(pd.read_csv(infile,sep='\t',encoding='utf-8_sig',usecols=col))
    
    #重排列
    df=pd.DataFrame(rgi_df,columns=col)
    #列标题重命名
    df.columns = ["比对基因ID","耐药基因ID","比对得分","置信度(%)","覆盖度(%)","耐药基因分类","作用机制"]  
    #排序
    df.sort_values(by=['比对得分','覆盖度(%)'],ascending=False,inplace=True)

    #读取深度file
    dep_df=pd.DataFrame(pd.read_csv(depfile,sep='\t',encoding='utf-8_sig'))
    dep_df.columns = ["比对基因ID","测序深度","基因覆盖度(%)"]

    #合并
    total_df=pd.merge(df,dep_df,how='left',on="比对基因ID")
    total_df['比对得分']=total_df['比对得分'].apply(lambda x: format(x, '.1f'))
    total_df['测序深度']=total_df['测序深度'].apply(lambda x: format(x, '.1f'))
    del total_df['基因覆盖度(%)']


    total_df.to_excel(f"{outdir}/detail_drug_resistance.xlsx",index=False,encoding='utf-8_sig')
    total_df.to_csv(f"{outdir}/detail_drug_resistance.txt",sep='\t',index=False,encoding='utf-8_sig')

else:
    exit(-1)