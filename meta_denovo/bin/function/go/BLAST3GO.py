#*coding=utf-8
#by chenweitian 2017
import optparse
import sys
import os
import gzip


usage="""
====================================================
***********BLAST3GO**********
Convert NCBI NR AccessionID to GO number
Author:chenweitian@genomics.cn
====================================================
"""
option = optparse.OptionParser(usage)
option.add_option('','--NR_tab',help='*result of BLAST NR,xml format',default='' )
option.add_option('','--ID_mapping',help='id mapping file',default='NA' )
option.add_option('','--GI_go',help='GI2GO',default='NA' )
option.add_option('','--Accession_go',help='accession2go,1th col=accessionId; 2nd=GO number',default='NA' )
option.add_option('','--output',help='*output',default='' )


(opts, args) = option.parse_args()
NR = opts.NR_tab
ID_MAP = opts.ID_mapping
OUT = opts.output
GI_GO = opts.GI_go
Accession_GO = opts.Accession_go

def main(NR,ID_MAP,OUT):
	Result = open(OUT,'w')
	GI2Gene = {}
	uniq = {}
	for line in open(NR,'r'):
		if line.startswith('Query_id'):continue
		line = line.strip().split('\t')
		if float(line[-3]) <= 1e-6:
			Gene = line[0]
			GI = line[1].split('|')[1]
			Key = "%s_%s"%(Gene,GI)
			if Key in uniq:
				continue
			else:
				uniq[Key] = 1
				GI2Gene.setdefault(GI,[]).append(Gene)
	Output_uniq = {}
	for line in open(ID_MAP,'r'):
		line = line.strip().split('\t')
		GI_MAP = line[4].split('; ')
		GO_MAP = line[7].split('; ')
		for Single_GI in GI_MAP:
			if Single_GI in GI2Gene:
				for Single_Gene in GI2Gene[Single_GI]:
					for Single_GO in GO_MAP:
						Key1 = "%s_%s"%(Single_GO,Single_GO)
						if Key1 in Output_uniq:
							continue
						else:
							Output_uniq[Key1] = 1
							Result.write("%s\t%s\n"%(Single_Gene,Single_GO))
	Result.close()

def main1(NR,GI_GO,OUT):
		Result = open(OUT,'w')
		GI2Gene = {}
		uniq = {}
		for line in open(NR,'r'):
			if line.startswith('Query_id'):continue
			line = line.strip().split('\t')
			if float(line[-3]) <= 1e-6:
				Gene = line[0]
				GI = line[1].split('|')[1]
				Key = "%s_%s"%(Gene,GI)
				if Key in uniq:
					continue
				else:
					uniq[Key] = 1
					GI2Gene.setdefault(GI,[]).append(Gene)

		Output_uniq = {}
		for line in open(GI_GO,'r'):
			line = line.strip().split('\t')
			GI_MAP = line[0]
			#print "%s\t%s"%(line[0],line[1])
			if len(line)<2:continue
			GO_MAP = line[1].split('; ')
		if GI_MAP in GI2Gene:
			for Single_Gene in GI2Gene[GI_MAP]:
				for Single_GO in GO_MAP:
					Result.write("%s\t%s\n"%(Single_Gene,Single_GO))
		Result.close()



def Accession2GO(NR_tab,DB,OUT):
	""" pass """
	Result = open(OUT,'w')
	Accession2Gene = {}
	uniqaccession = {}
	for line in open(NR,'r'):
		uniqaccession = {}
		if line.startswith('Query_id'):continue
		line = line.strip().split('\t')
		if float(line[-3]) <= 1e-6:
			Gene = line[0]

			accessionid0 = line[1]
			Accession2Gene.setdefault(accessionid0,set()).add(Gene)
			uniqaccession[accessionid0] = 1

			accessionid = line[-1].split('>')[1:]#split first ID which is accessionid0
			for accessionidi in accessionid:
				cluster_accessionID = accessionidi.split()[0]
				Accession2Gene.setdefault(cluster_accessionID,set()).add(Gene)
				

	with open(DB,'r') as ph:
		for line in ph:
			flag = 0 #cluster A & cluster B  = 1 
			line = line.strip().split('\t')
			if len(line)<2:continue
			AccessionID_col = line[0]
			GO_col = line[1]
			
			for AccessionID in AccessionID_col.split("; "):#NP_001193923.1; NP_001193924.1; NP_065149.2; NP_001193925.1; NP_870991.1
				if AccessionID in Accession2Gene:
					flag = 1
					for GeneID in Accession2Gene[AccessionID]:
						for GOID in GO_col.split("; "):
							#GO:0019888; GO:0005737; GO:0005856; GO:0016020; GO:0000159
							#print "%s\t%s"%(GeneID,GOID)
							Result.write("%s\t%s\n"%(GeneID,GOID))
				if flag == 1:
					break



if __name__ == '__main__':
	if len(sys.argv)<2:
		os.system("python %s -h"%(sys.argv[0]))
		sys.exit(1)

	else:
		if ID_MAP != "NA":
			main(NR,ID_MAP,OUT)
		elif GI_GO != "NA":
			main1(NR,GI_GO,OUT)
		elif Accession_GO != "NA":
			Accession2GO(NR,Accession_GO,OUT)

			
