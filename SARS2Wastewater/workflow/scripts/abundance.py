from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import matplotlib
matplotlib.use('Agg')


def getColor (var_name):
    rgbColors = plt.get_cmap('tab20b').colors + plt.get_cmap('tab20c').colors[0:16]
    colorCycle = [to_hex(color) for color in rgbColors]
    if var_name.lower() == 'other':
        return '#BBBBBB'
    else:
        color_idx = hash(var_name)%len(colorCycle)
        return colorCycle[color_idx]

def get_lineages_abundance_dict(demix):
    """谱系对应丰度的字典. {'XBB.1': 0.99, ...}"""
    dict_lngabdc = OrderedDict()
    with open(demix) as f:
        for line in f:
            lst = line.strip().split('\t')
            if lst[0] == 'lineages':
                lineages = lst[1].split(' ')
            elif lst[0] == 'abundances':
                abundances = lst[1].split(' ')
    for l,a in zip(lineages, abundances):
        dict_lngabdc[l] = float(a)
    return dict_lngabdc

def prepare4chart(dict_lngabdc):
    """获取饼图的输入数据"""
    thres_min_plot = 0.05
    names2plot, percentages2plot = [], []
    for lng in dict_lngabdc:
        if dict_lngabdc[lng] > thres_min_plot:
            names2plot.append(lng)
            percentages2plot.append(dict_lngabdc[lng])
    names2plot.append('Other')
    percentages2plot.append(1 - sum(percentages2plot))
    return names2plot, percentages2plot

def dray_pie_chart(names2plot, percentages2plot, outfilename) -> None:
    colors2plot = [getColor(name) for name in names2plot]
    explosionArray = np.full(len(percentages2plot), 0.07)
    plt.rcParams.update({'font.size': 12})
    plt.pie(percentages2plot, labels=names2plot, autopct='%1.1f%%', shadow=False,
            explode=explosionArray, colors=colors2plot)
    plt.axis('equal')
    plt.title('Abundance of variants')
    plt.savefig(outfilename, dpi=300)
    plt.close()

def generate_abundance_table(dict_lngabdc, abdctbl):
    """
    谱系分型丰度表格. 
    [230414]分型都列出来
    """
    with open(abdctbl, 'w') as g:
        g.write('分型\t丰度\n')
        for l in dict_lngabdc:
            g.write(f'{l}\t{dict_lngabdc[l]:.1%}\n')


if __name__ == '__main__':
    demix = snakemake.input[0]
    outchart = snakemake.output[0]
    abdctbl = snakemake.output[1]
    dict_lngabdc = get_lineages_abundance_dict(demix)
    names2plot, percentages2plot = prepare4chart(dict_lngabdc)
    dray_pie_chart(names2plot, percentages2plot, outchart)
    generate_abundance_table(dict_lngabdc, abdctbl)
    