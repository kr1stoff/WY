#!/usr/bin/python3
# @auther: stone
# @date: 2023-08-10
import sys
import pandas as pd
import os
import click


def assign_grade(row):
    if row['Con_percent'] >= 0.99 and row['Spe_per_same'] >= 99:
        if (row['Spe_same_min'] >= row['Spe_diff_max'] and row['Con_min_iden'] >= 99 and row['Con_min_cov'] >= 99):
            return 'S+'
        else:
            return 'S'
    elif row['Con_percent'] >= 0.97 and row['Spe_per_same'] >= 97:
        if (row['Spe_same_min'] >= row['Spe_diff_max'] and row['Con_min_iden'] >= 99 and row['Con_min_cov'] >= 99):
            return 'A+'
        else:
            return 'A'
    elif row['Con_percent'] >= 0.95 and row['Spe_per_same'] >= 95:
        if (row['Spe_same_min'] >= row['Spe_diff_max'] and row['Con_min_iden'] >= 99 and row['Con_min_cov'] >= 99):
            return 'B+'
        else:
            return 'B'
    elif row['Con_percent'] >= 0.92 and row['Spe_per_same'] >= 92:
        if (row['Spe_same_min'] >= row['Spe_diff_max'] and row['Con_min_iden'] >= 99 and row['Con_min_cov'] >= 99):
            return 'C+'
        else:
            return 'C'
    else:
        return 'D'


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input_result', required=True, type=click.Path(), help="输入文件,result格式")
@click.option('-o', '--output_result', required=True, type=click.Path(), help="输出,最后一列加上level")
def main(input_result, output_result):
    # 读取CSV文件
    result = os.path.abspath(input_result)
    df = pd.read_csv(result, sep='\t')

    # 在dfFrame中添加等级列
    df['level'] = df.apply(assign_grade, axis=1)
    df['product'] = df['Con_percent'] * df['Spe_per_same']
    sorted_df = df.sort_values(by=['product', 'Name'], ascending=[False, True])
    sorted_df = sorted_df.drop(columns=['product'])
    sorted_df.to_csv(output_result, sep='\t', index=False)


if __name__ == '__main__':
    main()
