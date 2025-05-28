#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import click
import sys
import logging
import pandas as pd
import math
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M', stream=sys.stdout)

# todo: 待修复BUG： 仅对rtab里面的列数进行计算，可能会存在某个基因组不满足cov/iden条件整个被过掉，少算了几个基因组，造成包容性偏大

# == 通用函数 =====================================================================


def get_file_content(file):
    out = os.path.dirname(os.path.abspath(file))
    title = str(os.path.basename(os.path.abspath(file))
                ).rsplit('.', maxsplit=1)[0]
    return out, title


def parse_file(table):
    col_number = pd.read_csv(table, sep='\t', encoding='utf-8').shape[1]
    df = pd.DataFrame(pd.read_csv(table, sep='\t', encoding='utf-8',
                      names=[i for i in range(1, col_number+1)]))
    logging.info(f"blast总共 {len(df.index)} 条记录")
    return df


def filter_by_column(df, col, val):
    filter_df = df[df[col] >= val]
    logging.info(f"过滤{val}后剩下 {len(filter_df.index)} 条记录")
    return filter_df


# === 插件功能 =====================================================================
@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def cli():
    """给定条件过滤blast结果"""
    pass


# == 覆盖度和置信度值过滤=====================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("blast_out", required=True, type=click.Path())
@click.option('-f', '--column', required=True, type=click.STRING, help="覆盖度和置信度列号，','分割")
@click.option('--cov', type=click.INT, default=95, help="覆盖度筛选，默认95%")
@click.option('--identity', type=click.INT, default=95, help="置信度筛选，默认95%")
@click.option('-o', '--out', required=True, type=click.Path(), help="输出文件")
def Filter(blast_out, column, cov, identity, out):
    """覆盖度和置信度过滤"""
    df = parse_file(blast_out)

    col_cov, col_iden = str(column).strip().split(',')[0:]
    df1 = filter_by_column(df, int(col_cov), int(cov))
    df2 = filter_by_column(df1, int(col_iden), int(identity))
    df2.to_csv(out, header=None, sep='\t', encoding='utf-8', index=False)
    logging.info(f"[Output]:  {out}")


# ===提取前N(条/%) blast记录=====================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("blast_out", required=True, type=click.Path())
@click.option('-n', '--top', type=click.FLOAT, default=10, help="覆盖度筛选，默认10，(小于1按百分比计算)")
@click.option('-o', '--out', required=True, type=click.Path(), help="输出文件")
def Top(blast_out, top, out):
    """提取前N(条/%) blast记录，第1列需为序列id"""
    top_df = pd.DataFrame()
    df = parse_file(blast_out)
    seqid_list = list(set(df.iloc[:, 0].values))
    logging.info(f"seq_id有 {len(seqid_list)} 条记录")

    for seq_id in seqid_list:
        tmp_df = df[df[1] == seq_id]
        if top >= 1:
            top = int(top)
            tmp_df = tmp_df.head(top)
        else:
            tmp_df = tmp_df.head(math.ceil(top*len(tmp_df)))
        top_df = pd.concat([top_df, tmp_df], axis=0,
                           join='outer', ignore_index=True)
    logging.info(f"top过滤后剩下 {len(top_df.index)} 条记录")

    top_df.to_csv(out, header=None, sep='\t', encoding='utf-8', index=False)
    logging.info(f"[Output]: {out}")


# ===统计比对基因数量 =====================================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("blast_out", required=True, type=click.Path())
@click.option('-f', '--column', required=True, type=click.STRING, help="query_id,ref_id,identity,cov列号，','分割")
@click.option('-o', '--out', required=True, type=click.Path(), help="输出文件")
def stat(blast_out, column, out):
    col_query, col_ref, col_iden, col_cov = str(column).strip().split(',')[0:]
    count_dic = {}
    cov_dic = {}
    iden_dic = {}
    with open(blast_out, 'r', encoding='utf-8') as f:
        for line in f:
            if not line:
                continue
            if line.startswith("qaccver"):
                continue

            query_id = line.strip().split('\t')[int(col_query)-1]
            ref_id = line.strip().split('\t')[int(col_ref)-1]
            identity = float(line.strip().split('\t')[int(col_iden)-1])
            cov = float(line.strip().split('\t')[int(col_cov)-1])

            count_dic.setdefault(query_id, {})
            count_dic[query_id].setdefault(ref_id, 0)
            count_dic[query_id][ref_id] += 1

            cov_dic.setdefault(query_id, [])
            cov_dic[query_id].append(cov)

            iden_dic.setdefault(query_id, [])
            iden_dic[query_id].append(identity)

    # Rtab矩阵
    Rtab_df = pd.DataFrame.from_dict(count_dic, orient='index')
    Rtab_df.fillna(0, inplace=True)       # 空白值0代替
    Rtab_df.insert(0, 'query_id', Rtab_df.index)
    Rtab_df.sort_values(by=['query_id'], ascending=True, inplace=True)
    Rtab_df.to_csv(out, index=False, sep='\t',
                   encoding='utf-8', float_format=int)
    logging.info(f"[Output]: {out}")

    # Rtab信息统计
    stat_df = pd.DataFrame()

    def count(col):
        tmp = len([i for i in list(col[1:]) if i >= 1])  # 计算 >1 的个数
        # tmp = int(sum(col[1:]))  # 累加
        return tmp

    # 计算 min/average
    def cal(col, dic, flag):
        query_id = col[int(col_query)-1]
        if flag == 'min':
            val = min(dic[query_id])
        elif flag == 'av':
            val = format(sum(dic[query_id])/len(dic[query_id]), '.2f')
        return val

    stat_df['count'] = Rtab_df.apply(count, axis=1)
    stat_df['total'] = len(Rtab_df.columns) - 1     # todo 需要减去标题？？
    stat_df['percent'] = (stat_df['count']/stat_df['total']
                          ).map(lambda x: format(x, '.2f'))
    stat_df['min_iden'] = Rtab_df.apply(cal, dic=iden_dic, flag='min', axis=1)
    stat_df['min_cov'] = Rtab_df.apply(cal, dic=cov_dic, flag='min', axis=1)
    stat_df['av_iden'] = Rtab_df.apply(cal, dic=iden_dic, flag='av', axis=1)
    stat_df['av_cov'] = Rtab_df.apply(cal, dic=cov_dic, flag='av', axis=1)

    stat_df.insert(0, 'query_id', stat_df.index)
    stat_df.sort_values(
        by=['percent', 'min_iden', 'min_cov'], ascending=False, inplace=True)
    stat_df.to_csv(f"{out}.stat", index=False, sep='\t', encoding='utf-8')
    logging.info(f"[Output]: {out}.stat")


# === 选项区 =====================================================================
cli.add_command(Filter)
cli.add_command(Top)
cli.add_command(stat)

if __name__ == '__main__':
    cli()
