
from encodings import utf_8_sig
import os
import sys  
import pandas as pd
import xlwt

if len(sys.argv) != 5:
    print('Usage: python %s <result_file> <des_file> <merge_bam.file> <out_dir>'%sys.argv[0])
    sys.exit(-1)

#传参
result_file=os.path.abspath(sys.argv[1])
des_file=os.path.abspath(sys.argv[2])
merge_bam=os.path.abspath(sys.argv[3])
out_dir=os.path.abspath(sys.argv[4])



def remake_result_file(result_file,des_file,merge_bam,out_dir):

    result_data=pd.read_csv(result_file, sep="\t")
    result_data=pd.DataFrame(result_data)
    result_data.columns=['比对氨基酸序列ID' ,'比对氨基酸序列长度' ,'毒力蛋白ID' ,'sgi', '参考蛋白质序列长度', '置信度(%)' ,'对齐氨基酸序列长度', '错配数' ,'gap数' ,'qstart', 'qend', 'sstart',' send', '期望值' ,'比对得分','覆盖度(%)', 'qcovhsp']
    
    result_data['覆盖度(%)']=result_data['覆盖度(%)'].apply(lambda x: int(x))
    result_data['比对得分']=result_data['比对得分'].astype(int)
    result_data['期望值']=result_data['期望值'].apply(lambda x: format(x, '.2e')) 
    print(f"total_len is {len(result_data)}")
    
    #去掉注释列 毒力蛋白ID 一样的
    result_data.drop_duplicates(subset=["毒力蛋白ID"],keep='first',inplace=True)
    print(f"\n\ndrop_duplicates total_len is {len(result_data)}")
    
    #条件筛选，覆盖度 >80
    cov=int(80)
    result_data=result_data.loc[(result_data["覆盖度(%)"] > cov)]
    print(f"cov>{cov} total_len has {len(result_data)}")
    

    #读入描述文件
    des_data=pd.read_csv(des_file,sep="\t")
    des_data=pd.DataFrame(des_data) 
    des_data.columns=['毒力蛋白ID', 'VF_Name', '毒力蛋白全称', '种类' ,'结构' ,'微生物体']

    #读入深度bed_bam_merge.depth文件
    merge_bam=pd.read_csv(merge_bam,sep="\t")
    merge_bam=pd.DataFrame(merge_bam) 
    merge_bam.columns=['比对氨基酸序列ID', '测序深度','seq_cov']
    merge_bam['测序深度']=merge_bam['测序深度'].apply(lambda x: format(x, '.1f'))


    #根据毒力蛋白ID，进行注释vlookup
    total_df=pd.merge(result_data,des_data,how='left',on="毒力蛋白ID")  
    total_df=pd.merge(total_df,merge_bam,how='left',on="比对氨基酸序列ID") 

    #去掉注释列 VF_FullName 一样的
    total_df.drop_duplicates(subset=["毒力蛋白全称"],keep='first',inplace=True)
    print(f"VF_FullName_dup total_len has {len(total_df)}")

    #删除含空值的行
    total_df.dropna(inplace=True)
    print(f"VF_FullName_Nan total_len has {len(total_df)}")

    #展示需要的列
    filter_df=total_df[['比对氨基酸序列ID','毒力蛋白ID','毒力蛋白全称','微生物体','覆盖度(%)','置信度(%)','参考蛋白质序列长度','比对氨基酸序列长度','对齐氨基酸序列长度','比对得分','测序深度']]

    total_df.to_csv(f'{out_dir}/detail_virulence_gene.txt',sep='\t',index=False,encoding='utf_8_sig')
    
    filter_df.to_csv(f'{out_dir}/show_virulence_gene.txt',sep='\t',index=False,encoding='utf_8_sig')
    filter_df.to_excel(f'{out_dir}/show_virulence_gene.xlsx',index=False,encoding='utf_8_sig')


#调用函数
if not os.path.exists(result_file) or not os.path.exists(merge_bam) or os.path.getsize(result_file) <int(5) or os.path.getsize(merge_bam) <int(5):
    
    a=open(f"{out_dir}/virulence_gene.txt",'w',encoding='utf_8_sig')
    a.write("比对氨基酸序列ID\t比对氨基酸序列长度\t毒力蛋白ID\tsgi\t参考蛋白质序列长度\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\t比对得分\t覆盖度(%)\tqcovhsp")   
    a.close()

    b=open(f"{out_dir}/show_virulence_gene.txt",'w',encoding='utf_8_sig')
    b.write("比对氨基酸序列ID\t毒力蛋白ID\t毒力蛋白全称\t微生物体\t覆盖度(%)\t置信度(%)\t参考蛋白质序列长度\t比对氨基酸序列长度\t对齐氨基酸序列长度\t比对得分\t测序深度")
    b.close()

    print(f"\nEmpty_file (merge_bam.file) or (virulence_gene.txt)")

else:
    remake_result_file(result_file,des_file,merge_bam,out_dir)


