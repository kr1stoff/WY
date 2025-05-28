#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/3/22 17:04
# @Last Modified by:   Ming
# @Last Modified time: 2022/3/22 17:04
import logging
import os

import click
import svgwrite
from Bio import SeqIO

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


#### Some Function
class Mutex(click.Option):
    """
    互斥选项参数
    """

    def __init__(self, *args, **kwargs):
        self.not_required_if: list = kwargs.pop("not_required_if")

        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs["help"] = (kwargs.get("help", "") + "\n[Option is mutually exclusive with " + ", ".join(
            self.not_required_if) + "]").strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt: bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '" + str(self.name) + "' is mutually exclusive with " + str(mutex_opt) + ".")
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)


class PrimerPosition(object):
    """
    Primer position plot
    """

    def __init__(self, record, f_out, width=None, height=None, top=10, bottom=10, left=10, right=10, gap=40):
        """
        Init the object

        :param record: The SeqRecord object of BioPython
        :param f_out: The out put file name
        :param width: The svg file width(default: None)
        :param height: The svg file height(default: None)
        :param top: The top margin(default: 10)
        :param bottom: The bottom margin(default: 10)
        :param left: The left margin(default: 10)
        :param right: The right margin(default: 10)
        :param gap: The gap between each part(default: 40)
        """
        self.out = os.path.abspath(f_out)
        self.record = record
        self.sequence = self.record.seq
        self.length = len(self.sequence)
        self.name = self.record.name
        self.top = top
        self.bottom = bottom
        self.left = left
        self.right = right
        self.gap = gap

        # 当前画笔的 x,y
        self.x = self.left
        self.y = self.top
        self._unit()

        # 使用等宽字体monospace确保所有字符宽度一直
        # 5号字体宽度为2.5，长度为
        self.font = "monospace"
        self.fontsize = 4
        self.fontwidth = 2
        self.fontheight = 4

        self.width = width if width else self.length * self.fontwidth + left + right
        self.height = height if height else 500

        self.svg = svgwrite.Drawing(f_out,
                                    size=(self.width, self.height),
                                    debug=True)

        self._add_chromosome()

    def _add_chromosome(self):
        """
        Draw the chromosome line
        """
        # 染色体高度
        width = self.length * self.fontwidth
        # 染色体
        self.svg.add(self.svg.rect((self.x, self.y),
                                   (width, self.fontheight),
                                   stroke="grey",
                                   fill='none'))
        # 碱基
        # 碱基字体大小
        self.svg.add(self.svg.text(self.sequence,
                                   insert=(self.left, self.y + self.fontheight),
                                   style=f"text-anchor:start;baseline-shift:10%;font-size:{self.fontsize}px; font-family:{self.font}"))
        self.y += self.fontheight

        # 主刻度
        tickheight = 2.5
        for i in range(0, self.length // self.unit + 1):
            self.svg.add(self.svg.line(start=(self.left + self.unit * i * self.fontwidth, self.y),
                                       end=(self.left + self.unit * i * self.fontwidth, self.y + tickheight),
                                       stroke="black"))
            if self.unit == 1000:
                text = f"{i}K"
            elif self.unit == 1000000:
                text = f"{i}M"
            else:
                text = i * self.unit
            self.svg.add(self.svg.text(text,
                                       insert=(self.left + self.unit * i * self.fontwidth, self.y + tickheight),
                                       style=f"text-anchor:middle;baseline-shift:-80%;font-size:{self.fontsize}px; font-family:{self.font}"))
        # 副刻度
        for i in range(0, self.length // (self.unit // 10) + 1):
            self.svg.add(self.svg.line(start=(self.left + self.unit // 10 * i * self.fontwidth, self.y),
                                       end=(self.left + self.unit // 10 * i * self.fontwidth, self.y + tickheight / 2),
                                       stroke="grey"))

        self.y += tickheight
        self.y += self.gap

    def plot_snp(self, tag, pos, y=None):
        """
        Plot the snp site
        """
        y = y if y else self.y
        self.svg.add(self.svg.text(tag,
                                   insert=(self.left + pos[0] * self.fontwidth, y),
                                   style=f"text-anchor:middle;baseline-shift:-33%;font-size:{self.fontsize}px; font-family:{self.font}"))

    def plot_background(self, color, opacity=0.5):
        """
        Plot a track background
        """
        self.svg.add(self.svg.rect((self.x, self.y - self.fontheight / 2),
                                   (self.length * self.fontwidth, self.fontheight),
                                   stroke="none",
                                   opacity=opacity,
                                   fill=color))

    def plot_gene(self):
        """
        Plot the gene element
        """
        colors = ["#A1C5DE", "#F6CEA8", "#EFB8BA", "#DBF9DF", "#D0BEE2", "#B7B7B8"]
        height = 10

        flag = 0
        for i in self.record.features:
            if i.type == 'gene':
                color = colors[flag % len(colors)]
                gene_name = i.qualifiers["gene"][0]
                pos = (int(i.location.start), int(i.location.end))

                self.svg.add(self.svg.rect((self.left + pos[0] * self.fontwidth, self.y),
                                           ((pos[1] - pos[0]) * self.fontwidth, height),
                                           stroke="grey",
                                           fill=color))
                self.svg.add(self.svg.text(gene_name,
                                           insert=((pos[0] + pos[1]) * self.fontwidth / 2, self.y - 5),
                                           stroke="black",
                                           style=f"text-anchor:middle;baseline-shift:-33%;font-size:10px; font-family:{self.font}"))
                flag += 1
        self.y += (height + self.gap)

    def plot_snp(self, pos, base, color="red", y=None):
        """
        Plot the snp info
        """
        y = y if y else self.y

        self.svg.add(self.svg.text(base,
                                   insert=(self.left + pos * self.fontwidth, y),
                                   style=f"text-anchor:start;baseline-shift:-33%;font-size:{self.fontsize}px; font-family:{self.font}"))

    def plot_primer(self, primers):
        """
        Plot the primer element

        :param primers: The primer info list
        """

        def is_overlap(primer1, primer2):
            """
            Check whether primer1 and primer2 has overlap
            """
            if primer1 is None or primer2 is None:
                return False

            p1_start, p1_end = primer1[1], primer1[4]
            p2_start, p2_end = primer2[1], primer2[4]
            if (p1_start <= p2_start <= p1_end) or (p1_start <= p2_end <= p1_end):
                return True
            else:
                return False

        # 有重叠的引物的偏差
        ybias = 10
        # 产物高度
        height = 10
        last_primer = None
        last_y = self.y
        for primer in primers:
            if is_overlap(last_primer, primer):
                y = last_y + ybias
            else:
                y = self.y

            # left
            self.svg.add(self.svg.line(start=(self.left + primer[1] * self.fontwidth, y),
                                       end=(self.left + primer[2] * self.fontwidth, y),
                                       stroke="#c7001f"))
            self.svg.add(self.svg.text(primer[5],
                                       insert=(self.left + primer[1] * self.fontwidth, y),
                                       style=f"text-anchor:start;baseline-shift:-33%;font-size:{self.fontsize}px; font-family:{self.font}"))
            # right
            self.svg.add(self.svg.line(start=(self.left + primer[3] * self.fontwidth, y),
                                       end=(self.left + primer[4] * self.fontwidth, y),
                                       stroke="#c7001f"))
            self.svg.add(self.svg.text(primer[6],
                                       insert=(self.left + primer[3] * self.fontwidth, y),
                                       style=f"text-anchor:start;baseline-shift:-33%;font-size:{self.fontsize}px; font-family:{self.font}"))

            # product
            self.svg.add(self.svg.line(start=(self.left + primer[2] * self.fontwidth, y),
                                       end=(self.left + primer[3] * self.fontwidth, y),
                                       stroke="green",
                                       stroke_width=4))
            last_y = y
            last_primer = primer

        self.y += height
        self.y += self.gap

    def draw(self):
        """
        Generate the out put svg file
        """
        self.svg.save()

    def _unit(self):
        """
        Judge the unit to use for tick
        """
        self.unit = 100
        if self.length > 10000:
            self.unit = 1000
        if self.length > 10000000:
            self.unit = 1000000


##################
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--genbank",
              cls=Mutex,
              not_required_if=["fasta"],
              type=click.Path(),
              help="The input genbank format genome file")
@click.option("--fasta",
              cls=Mutex,
              not_required_if=["genbank"],
              type=click.Path(),
              help="The input fasta format genome file")
@click.option('--primer',
              required=True,
              multiple=True,
              type=click.Path(),
              help="The primer info file")
@click.option('--vcf',
              multiple=True,
              type=click.Path(),
              help="The snp vcf info file(现在只是一个突变位点的tab文件，并非vcf格式文件)")
@click.option('--height',
              default=500,
              show_default=True,
              help="The svg height")
@click.option('--gap',
              default=40,
              show_default=True,
              help="The gap of tracks")
@click.option('--out',
              default="./",
              show_default=True,
              type=click.Path(),
              help="The out put dir")
def main(genbank, fasta, primer, vcf, height, gap, out):
    """
    Visualization of primers info on the genome

    引物格式（无头，最后一列可选）
    染色体\tF起点\tF终点\tR起点\tR终点\tF序列\tR序列\t名称
    """
    if len(vcf) > 0:
        logging.info(f"Parse the mutation info")
        info_snp = {}
        for f_vcf in vcf:
            with open(f_vcf, 'r') as IN:
                for line in IN:
                    arr = line.strip().split('\t')
                    info_snp[arr[0]] = [[i[0], int(i[1:-1]), i[-1]] for i in arr[1].strip().split(';')]

    logging.info(f"Start to draw")
    if genbank:
        fformat = "genbank"
        f_genome = os.path.abspath(genbank)
    else:
        fformat = "fasta"
        f_genome = os.path.abspath(fasta)

    for record in SeqIO.parse(f_genome, fformat):
        name = record.name
        fout = os.path.join(out, f"{name}.svg")

        res = PrimerPosition(record, fout, height=height, gap=int(gap))
        if fformat == "genbank":
            res.plot_gene()

        # TODO: snp的绘制待确定输入的格式后集成到类中
        # 现有输入格式
        # Alpha   A29782T;C17678T;A29T
        # Beta    C15714T;C9142T;G23012A
        colors = ["#A1C5DE", "#F6CEA8", "#EFB8BA", "#DBF9DF", "#D0BEE2", "#B7B7B8"]
        flag = 0
        for k, v in info_snp.items():
            res.plot_background(colors[flag])
            for pos in v:
                res.plot_snp(pos[1], pos[2], y=res.y)
            res.y += 5
            flag += 1
            flag = flag % len(colors)
        res.y += int(gap) - 5

        # 引物
        for f_primer in primer:
            info_primer = []
            with open(f_primer, 'r') as IN:
                for line in IN:
                    if line.startswith("#"):
                        pass
                    else:
                        arr = line.strip().split('\t')
                        arr[1:5] = [int(i) for i in arr[1:5]]
                        info_primer.append(arr)
            info_primer.sort(key=lambda x: x[1])
            res.plot_primer(info_primer)

        res.draw()


if __name__ == "__main__":
    main()
