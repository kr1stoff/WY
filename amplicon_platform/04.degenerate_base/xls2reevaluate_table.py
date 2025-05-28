import sys

xls = sys.argv[1]

with open(xls, 'r') as IN:
    for line in IN:
        if line.startswith('species'):
            continue
        else:
            arr = line.split('\t')
            name = arr[5]
            Left_Seq = arr[2]
            Right_Seq = arr[3]
            Probe_Seq = arr[4]
            Left_Seq_out = arr[29]
            Right_Seq_out = arr[30]
            Taxid = arr[46]
            name2 = name.replace(' ', '_')
            print(name2, Left_Seq, Right_Seq, Probe_Seq, Left_Seq_out,
                  Right_Seq_out, Taxid, Taxid, sep='\t')
