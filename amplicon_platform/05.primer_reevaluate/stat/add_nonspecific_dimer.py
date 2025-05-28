import click


def read_table(table):
    dict_table = {}
    with open(table, 'r') as IN:
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            name = arr[0]
            content = '\t'.join(arr[1:])
            dict_table[name] = content
    return dict_table


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-s', '--stat', required=True, type=click.Path(), help="特异性,包容性结果表格")
@click.option('-n', '--nonspecific', required=True, type=click.Path(), help="非特异性扩增表格")
@click.option('-do', '--dimer_in', required=True, type=click.Path(), help="内引物二聚体结果")
@click.option('-do', '--dimer_out', required=True, type=click.Path(), help="外引物二聚体结果")
@click.option('-o', '--out', default='./primer.merge.txt', type=click.Path(), help="输出")
def main(stat, nonspecific, dimer_in, dimer_out, out):
    dict_stat = read_table(stat)
    dict_nonspecific = read_table(nonspecific)
    dict_dimer_in = read_table(dimer_in)
    dict_dimer_out = read_table(dimer_out)
    with open(stat, 'r') as IN, open(out, 'w') as OUT:
        for line in IN:
            line = line.strip()
            arr = line.split('\t')
            if line.startswith('ID'):
                print(
                    f'{line}\tNonspecific\tInprimer_dimer_F\tInprimer_dimer_R\tOutprimer_dimer_F\tOutprimer_dimer_R', file=OUT)
            else:
                name = arr[0]
                dict_stat.setdefault(name, '-\t')
                dict_nonspecific.setdefault(name, '-')
                dict_dimer_in.setdefault(name, '-\t-')
                dict_dimer_out.setdefault(name, '-\t-')
                print(
                    f'{name}\t{dict_stat[name]}\t{dict_nonspecific[name]}\t{dict_dimer_in[name]}\t{dict_dimer_out[name]}', file=OUT)


if __name__ == '__main__':
    main()
