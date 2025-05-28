import click


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input', required=True, type=click.Path(), help="输入")
@click.option('-t', '--table', required=True, type=click.Path(), help="带有ID的表格")
@click.option('-s', '--suffixes', required=True, type=click.Path(), help="ID的后缀")
@click.option('-o', '--output', default='./out.txt', type=click.Path(), help="输出")
def main(input, table, suffixes, output):
    dict = {}
    with open(input, 'r') as IN:
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            dimer_pair = arr[1]
            score = arr[2]
            length = arr[9]
            mismatch = arr[10]

            pairs = dimer_pair.split('::')

            if (int(score) >= 6):
                dict.setdefault(pairs[0], [])
                dict.setdefault(pairs[1], [])
                dict[pairs[0]].append(f'{pairs[1]}|{length}-{mismatch}')
                dict[pairs[1]].append(f'{pairs[0]}|{length}-{mismatch}')

    with open(table, 'r') as P, open(output, 'w') as OUT:
        for line in P:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            f = name + "__F__" + suffixes
            r = name + "__R__" + suffixes
            if (f in dict.keys() and r in dict.keys()):
                f_out = ",".join(dict[f])
                r_out = ",".join(dict[r])
                print(f'{name}\t{f_out}\t{r_out}', file=OUT)
            elif (f in dict.keys()):
                f_out = ",".join(dict[f])
                print(f'{name}\t{f_out}\t-', file=OUT)
            elif (r in dict.keys()):
                r_out = ",".join(dict[r])
                print(f'{name}\t-\t{r_out}', file=OUT)
            else:
                print(f'{name}\t-\t-', file=OUT)


if __name__ == '__main__':
    main()
