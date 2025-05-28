import sys

taxid_list = sys.argv[1]
merge = sys.argv[2]

background_list = []
with open(taxid_list, 'r') as IN:
    for line in IN:
        arr = line.split('\t')
        background_list.append(arr[0])

with open(merge, 'r') as IN:
    taxids = []
    for line in IN:
        arr = line.strip().split('\t')
        in_same_diff = arr[22]
        out_same_diff = arr[23]
        probe_same_diff = arr[24]
        name = arr[0]
        taxid = arr[2]
        if (in_same_diff == "Diff" and out_same_diff == "Diff" and probe_same_diff == "Diff" and taxid in background_list):
            taxids.append(taxid)
    unique_set = set(taxids)
    taxid_list = ','.join(unique_set)
    print(f'{name}\t{taxid_list}')
