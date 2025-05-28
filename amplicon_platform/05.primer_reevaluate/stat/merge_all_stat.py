#!/usr/bin/python3
# @auther: stone
# @date: 2023-07-26
import sys

input = sys.argv[1]
specificity_in = sys.argv[2]
specificity_p = sys.argv[3]
specificity_out = sys.argv[4]
specificity_merge = sys.argv[5]

inclusiveness_in = sys.argv[6]
inclusiveness_p = sys.argv[7]
inclusiveness_out = sys.argv[8]
inclusiveness_merge = sys.argv[9]


def get_result(file, col):
    result = {}
    with open(file, 'r') as IN:
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            per = arr[(int(col) - 1)]
            result[name] = per
    return result


def main(input, specificity_in, specificity_p, specificity_out, specificity_merge, inclusiveness_in, inclusiveness_p, inclusiveness_out, inclusiveness_merge):
    dic_specificity_in = get_result(specificity_in, 2)
    dic_specificity_p = get_result(specificity_p, 2)
    dic_specificity_out = get_result(specificity_out, 2)
    dic_specificity_merge = get_result(specificity_merge, 2)

    dic_inclusiveness_in = get_result(inclusiveness_in, 2)
    dic_inclusiveness_p = get_result(inclusiveness_p, 2)
    dic_inclusiveness_out = get_result(inclusiveness_out, 2)
    dic_inclusiveness_merge = get_result(inclusiveness_merge, 4)
    with open(input, 'r') as IN, open('all.primer.merge.stat', 'w') as OUT:
        print(*['ID', 'In_specificity', 'Probe_specificity', 'Out_specificity', 'In_inclusiveness', 'Probe_inclusiveness',
              'Out_inclusiveness', 'Merge_specificity', 'Merge_inclusiveness'], sep='\t', file=OUT)
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            dic_specificity_in.setdefault(name, 0)
            dic_specificity_p.setdefault(name, 0)
            dic_specificity_out.setdefault(name, 0)
            dic_inclusiveness_in.setdefault(name, 0)
            dic_inclusiveness_p.setdefault(name, 0)
            dic_inclusiveness_out.setdefault(name, 0)
            dic_specificity_merge.setdefault(name, 0)
            dic_inclusiveness_merge.setdefault(name, 0)
            print(*[name, dic_specificity_in[name], dic_specificity_p[name], dic_specificity_out[name], dic_inclusiveness_in[name],
                  dic_inclusiveness_p[name], dic_inclusiveness_out[name], dic_specificity_merge[name], dic_inclusiveness_merge[name]], sep='\t', file=OUT)


if __name__ == '__main__':
    main(input, specificity_in, specificity_p, specificity_out, specificity_merge,
         inclusiveness_in, inclusiveness_p, inclusiveness_out, inclusiveness_merge)
