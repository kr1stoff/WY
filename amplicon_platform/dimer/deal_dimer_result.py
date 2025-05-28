import sys
import re

dimer = sys.argv[1]
primer = sys.argv[2]

dict = {}
with open(dimer, 'r') as IN:
    for line in IN:
        arr = line.strip().split('\t')
        dimer_pair = arr[1]
        score = arr[2]
        length = arr[9]
        mismatch = arr[10]
        flag = arr[11]

        pairs = dimer_pair.split('::')

        if (flag == "Y" and int(length) >= 7 and int(mismatch) <= 1 and int(score) >= 6):
            dict.setdefault(pairs[0], [])
            dict.setdefault(pairs[1], [])
            if re.search(r'\bS[01]', pairs[1]):
                dict[pairs[0]].append(f'{pairs[1]}|{length}-{mismatch}')
            if re.search(r'\bS[01]', pairs[0]):
                dict[pairs[1]].append(f'{pairs[0]}|{length}-{mismatch}')

# print (dict)

with open(primer, 'r') as P:
    for line in P:
        line = line.strip()
        f = line + "_F"
        r = line + "_R"
        if (f in dict.keys() and dict[f]):
            print(f'{line}\t{",".join(dict[f])}')
        elif (r in dict.keys() and dict[r]):
            print(f'{line}\t{",".join(dict[r])}')
        else:
            print(f'{line}\t-')
