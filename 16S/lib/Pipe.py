#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/1/18 15:49
# @Last Modified by:   Ming
# @Last Modified time: 2022/1/18 15:49
import logging
import os
import os.path
import sys

import yaml

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')
LOG = logging.getLogger(__name__)

#### My own lib
BIN = os.path.dirname(os.path.abspath(__file__))
d_bin = os.path.join(BIN, "../bin")
d_config = os.path.join(BIN, "../config")
f_yaml = os.path.join(d_config, "software.yml")
f_env = os.path.join(d_config, "env.yml")
f_database = os.path.join(d_config, "database.yml")
config = yaml.safe_load(open(f_yaml, 'r').read())
# PATH_EasyPipe = config["EasyPipe"]
# sys.path.append(PATH_EasyPipe)
from EasyPipe import Pipe


class Meta(Pipe):
    """
    16S Pipeline for EARTH
    """

    def __init__(self, f_config, d_out):
        """
        Init the Pipe
        :param f_config: The config for 16S
        :param d_out: The out put dir for 16S
        """
        self.f_config = os.path.abspath(f_config)
        self.out = os.path.abspath(d_out)

        self.config = yaml.safe_load(open(self.f_config, 'r').read())
        self.name = self.config["project"]
        super(Meta, self).__init__(self.name,
                                   self.out,
                                   bin_folder=d_bin,
                                   config_folder=d_config)
        self.rawdata = self.config["rawdata"]

        self.env = yaml.safe_load(open(f_env, 'r').read())
        self.database = yaml.safe_load((open(f_database, 'r').read()))

        # 一些流程种经常用到的变量
        self.samples = self.config["samples"]
        self.groups = self.config[
            "group_order"] if "group_order" in self.config else None
        if self.groups:
            self.sample2group = {}
            for group in self.groups:
                for sample in self.config["groups"][group]:
                    self.sample2group[sample] = group
            self.group_diff = self.config["group_diff"] if "group_diff" in self.config else None
        else:
            self.sample2group = None
            self.group_diff = None

        # 常用软件
        self.qiime = self.software["qiime2"]
        self.python = self.software["python"]
        self.perl = self.software["perl"]
        self.rscript = self.software["Rscript"]
        self.biom = self.software["biom"]
        self.generalplot = f"{self.perl} {self.software['GeneralPlot']}"

        # 需要剪切掉的长度
        self.trim = {("16S", "V3-V4"): (17, 20),
                     ("16S", "V4"): (19, 20),
                     ("ITS", "ITS2"): (18, 20),
                     ("18S", "V4"): (19, 20)}

    def prepare(self):
        """
        Prepare files for qiime2
        """
        dir_out = os.path.join(self.out, "Prepare")
        job = self.child("Prepare", dir_out, _type="Job")

        # 准备 manifest.tsv 文件
        f_manifest = os.path.join(dir_out, "manifest.tsv")
        self.generate_manifest(self.config["samples"],
                               self.config["rawdata"],
                               f_manifest,
                               data_type=self.config["data_type"])
        # 准备 metadata.tsv 文件
        if self.sample2group:
            f_metadata = os.path.join(dir_out, "metadata.tsv")
            self.generate_metadata(f_metadata)

        # 准备group.qiime.tsv 文件
        if self.sample2group:
            f_qiimegroup = os.path.join(dir_out, "group.qiime.tsv")
            self.generate_qiimegroup(f_qiimegroup)

    def qc(self):
        """
        The QC step
        """
        dir_out = os.path.join(self.out, "QC")
        job = self.child("QC", dir_out)

        d_fastqc = os.path.join(dir_out, "FastQC")
        j_fastqc = job.child("fastqc", d_fastqc, run_type="multi_run",
                             _type="Job")
        j_fastqc.set_maxjob(self.config["parallel"])
        fastqc = self.software["fastqc"]
        for sample in self.samples:
            f1 = os.path.join(self.config["rawdata"], f"{sample}_1.fq.gz")
            f2 = os.path.join(self.config["rawdata"], f"{sample}_2.fq.gz")
            cmd = f"{fastqc} --quiet --extract -t 2 -o {d_fastqc} {f1} {f2} "
            j_fastqc.add_command(cmd)

    def qiime2(self):
        """
        QIIME2 for 16S
        可视化文件可倒入 https://view.qiime2.org 查看
        """
        dir_out = os.path.join(self.out, "QIIME2")
        job = self.child("QIIME2", dir_out)
        itype = {"PE": 'SampleData[PairedEndSequencesWithQuality]',
                 "SE": 'SampleData[SequencesWithQuality]'}
        iformat = {("PE", 33): "PairedEndFastqManifestPhred33V2",
                   ("PE", 64): "PairedEndFastqManifestPhred64V2",
                   ("SE", 33): "SingleEndFastqManifestPhred33V2",
                   ("SE", 64): "SingleEndFastqManifestPhred64V2"}

        # 向qiime2中导入数据
        f_manifest = os.path.join(self.out, "Prepare/manifest.tsv")
        d_import = os.path.join(dir_out, "Data")
        j_import = job.child("Import", d_import, _type="Job")
        data_type = self.config["data_type"]
        data_quality = self.config["data_quality"]
        j_import.add_command(self.env_start(self.env["qiime2"]))
        cmd = f"{self.qiime} tools import --type {itype[data_type]} --input-path {f_manifest} --output-path {d_import}/data.qza --input-format {iformat[(data_type, data_quality)]}"
        j_import.add_command(cmd)

        # Denois
        d_denois = os.path.join(dir_out, "Denois")
        j_denois = job.child("Denois", d_denois, _type="Job")
        j_denois.add_command(self.env_start(self.env["qiime2"]))
        # TODO: 后续需要更改获取序列的方式，是否需要从 qiime2 中独立出来，前端需要添加 qiime2 参数
        cmd = f"""
{self.qiime} dada2 denoise-paired --i-demultiplexed-seqs {d_import}/data.qza --p-trunc-len-f 0 --p-trunc-len-r 0 --o-table {d_denois}/table.qza --o-representative-sequences {d_denois}/rep-seqs.qza --o-denoising-stats {d_denois}/denoising-stats.qza --p-n-threads {self.config['threads']}
# {self.qiime} feature-table filter-features --i-table {d_denois}/table.qza --p-min-frequency 3 --o-filtered-table {d_denois}/table.filter.qza
{self.qiime} tools export --input-path {d_denois}/denoising-stats.qza --output-path {d_denois}/stat
sed -i 2d {d_denois}/stat/stats.tsv
{self.software["dos2unix"]}  {d_denois}/stat/stats.tsv
{self.qiime} tools export --input-path {d_denois}/table.qza --output-path {d_denois}/feature
{self.biom} convert -i {d_denois}/feature/feature-table.biom -o {d_denois}/feature/feature-table.tsv --to-tsv
{self.qiime} tools export --input-path {d_denois}/rep-seqs.qza --output-path {d_denois}/fasta
{self.python} {self.bin}/FeatureRename.py -f {d_denois}/feature/feature-table.tsv -s {d_denois}/fasta/dna-sequences.fasta -o {d_denois}/rename
{self.qiime} tools import --input-path {d_denois}/rename/dna-sequences.fasta --output-path {d_denois}/rename/rep-seqs.qza --type 'FeatureData[Sequence]'
{self.biom} convert -i {d_denois}/rename/feature-table.tsv -o {d_denois}/rename/feature-table.biom --table-type 'OTU table' --to-hdf5
{self.qiime} tools import --input-path {d_denois}/rename/feature-table.biom --output-path {d_denois}/rename/feature.qza --type 'FeatureTable[Frequency]' --input-format BIOMV210Format
sed 1d {d_denois}/rename/feature-table.tsv | sed 1s',#OTU ,,'g > {d_denois}/rename/feature.table.xls
{self.env_end()}
{self.python} {self.bin}/DenoisStat.py -f {d_denois}/stat/stats.tsv -o {d_denois}/stat/plot.data
{self.generalplot} bar -file {d_denois}/stat/plot.data -outprefix {d_denois}/stat/denois.fill -barfill T
{self.generalplot} bar -file {d_denois}/stat/plot.data -outprefix {d_denois}/stat/denois.stack -barfill F
"""
        j_denois.add_command(cmd)

        # phylogeny
        d_tree = os.path.join(dir_out, "Tree")
        j_tree = job.child("Tree", d_tree, _type="Job")
        j_tree.add_command(self.env_start(self.env["qiime2"]))
        cmd = f"""
{self.qiime} alignment mafft --i-sequences {d_denois}/rename/rep-seqs.qza --o-alignment {d_tree}/aligned-rep-seqs.qza
{self.qiime} alignment mask --i-alignment {d_tree}/aligned-rep-seqs.qza --o-masked-alignment {d_tree}/masked-aligned-rep-seqs.qza
{self.qiime} phylogeny fasttree --i-alignment {d_tree}/masked-aligned-rep-seqs.qza --o-tree {d_tree}/unrooted-tree.qza
{self.qiime} phylogeny midpoint-root --i-tree {d_tree}/unrooted-tree.qza --o-rooted-tree {d_tree}/rooted-tree.qza
{self.qiime} tools export --input-path {d_tree}/rooted-tree.qza --output-path {d_tree}/rooted
"""
        j_tree.add_command(cmd)

    def taxon(self):
        """
        Taxon analysis
        """
        dir_out = os.path.join(self.out, "Taxon")
        job = self.child("Taxon", dir_out)

        d_prepare = os.path.join(self.out, "Prepare")
        d_qiime2_result = os.path.join(self.out, "QIIME2/Denois/rename")
        align_database = self.config["database"]
        j_taxon = job.child("taxon", dir_out, _type="Job")
        j_taxon.add_command(self.env_start(self.env["qiime2"]))
        cmd = f"""
{self.qiime} feature-classifier classify-sklearn \\
    --i-classifier {self.database[align_database]} \\
    --i-reads {d_qiime2_result}/rep-seqs.qza \\
    --o-classification {dir_out}/taxonomy.qza
{self.qiime} metadata tabulate \\
    --m-input-file {dir_out}/taxonomy.qza \\
    --o-visualization {dir_out}/taxonomy.qzv
{self.qiime} tools export \\
    --input-path {dir_out}/taxonomy.qza \\
    --output-path {dir_out}/taxonomy
{self.python} {self.bin}/FormatLineage.py -i {dir_out}/taxonomy/taxonomy.tsv -n Taxon --new_sep ';' -o {dir_out}/taxonomy/taxonomy.format.tsv
{self.qiime} taxa barplot \\
    --i-table {d_qiime2_result}/feature.qza \\
    --i-taxonomy {dir_out}/taxonomy.qza \\
    --m-metadata-file {d_prepare}/metadata.tsv \\
    --o-visualization {dir_out}/taxa-bar-plots.qzv
{self.qiime} taxa collapse \\
    --i-table {d_qiime2_result}/feature.qza \\
    --i-taxonomy {dir_out}/taxonomy.qza \\
    --p-level 7 \\
    --o-collapsed-table {dir_out}/feature.taxonomy.qza
{self.qiime} feature-table relative-frequency \\
    --i-table {dir_out}/feature.taxonomy.qza \\
    --o-relative-frequency-table {dir_out}/feature.taxonomy.relative.qza
{self.qiime} tools export \\
    --input-path {dir_out}/feature.taxonomy.relative.qza \\
    --output-path {dir_out}/feature.taxonomy.relative
{self.biom} convert \\
    -i {dir_out}/feature.taxonomy.relative/feature-table.biom \\
    -o {dir_out}/feature.taxonomy.relative/feature-table.xls \\
    --to-tsv
sed -i 1d {dir_out}/feature.taxonomy.relative/feature-table.xls
sed -i "1s,#OTU ID,Lineage,"g {dir_out}/feature.taxonomy.relative/feature-table.xls
{self.python} {self.bin}/FormatLineage.py -i {dir_out}/feature.taxonomy.relative/feature-table.xls -n Lineage -o {dir_out}/feature.taxonomy.relative/feature-table.format.xls
"""
        j_taxon.add_command(cmd)

        # Taxon Stat
        d_stat = os.path.join(dir_out, "TaxonStat")
        j_stat = job.child("taxon_stat", d_stat, _type="Job")
        cmd = f"""
{self.perl} {self.bin}/make_qiime_table.pl {d_qiime2_result}/feature.table.xls {dir_out}/taxonomy/taxonomy.format.tsv > {d_stat}/all.feature.tab
{self.rscript} {self.bin}/Otu.recast.R -i {d_stat}/all.feature.tab -o {d_stat}/TaxonProfile
{self.generalplot} bar -file {d_stat}/TaxonProfile/all.taxonomy.stat.xls -outprefix {d_stat}/TaxonProfile/taxonomy.fill -barfill T
{self.generalplot} bar -file {d_stat}/TaxonProfile/all.taxonomy.stat.xls -outprefix {d_stat}/TaxonProfile/taxonomy.stack -barfill F
"""
        j_stat.add_command(cmd)
        levels = ["Domain", "Phylum", "Order", "Family", "Genus", "Species"]
        for i in levels:
            f_name = os.path.join(d_stat,
                                  f"TaxonProfile/profiling.{i}.top10.xls")
            cmd = f"{self.perl} {self.bin}/dfplot.pl {f_name} t {self.software['GeneralPlot']} bar -barstack T -outprefix {d_stat}/TaxonProfile/{i}.stack"
            j_stat.add_command(cmd)
        for i in levels[1:]:
            f_name = os.path.join(d_stat,
                                  f"TaxonProfile/profiling.{i}.min0.1.xls")
            cmd = f"{self.generalplot} cluster -file {f_name} -outprefix {d_stat}/TaxonProfile/{i}.heatmap -clustrow T -clustcol F"
            j_stat.add_command(cmd)

    def diversity(self):
        """
        Diversity analysis
        """
        # TODO: 把Alpha移动到单独的方法
        dir_out = os.path.join(self.out, "Diversity")
        job = self.child("Diversity", dir_out, run_type="multi_run")
        job.set_maxjob(3)

        d_prepare = os.path.join(self.out, "Prepare")
        d_qiime2 = os.path.join(self.out, "QIIME2")
        d_tree = os.path.join(d_qiime2, "Tree")

        #         # 计算核心多样性
        #         d_core = os.path.join(dir_out, "Core")
        #         j_core = job.child("Core", d_core, _type="Job")
        #         cmd = f"{self.python} {self.bin}/GetMinDepth.py -f {d_qiime2}/Denois/rename/feature.table.xls -o {dir_out}/min.depth.txt"
        #         j_core.add_command(cmd)
        #         j_core.add_command(self.env_start(self.env["qiime2"]))
        #         cmd = f"""
        # if [ -d {d_core} ]; then
        #     rm -rf {d_core}
        # fi
        # cat {dir_out}/min.depth.txt | \\
        # {self.qiime} diversity core-metrics-phylogenetic \\
        #             --i-phylogeny {d_tree}/rooted-tree.qza \\
        #             --i-table {d_qiime2}/Denois/rename/feature.qza \\
        #             --p-sampling-depth - \\
        #             --m-metadata-file {d_prepare}/metadata.tsv \\
        #             --output-dir {d_core} \\
        #             --p-n-jobs-or-threads {self.config["threads"]}
        # {self.qiime} diversity alpha-correlation \\
        #             --i-alpha-diversity {d_core}/observed_features_vector.qza \\
        #             --m-metadata-file {d_prepare}/metadata.tsv \\
        #             --o-visualization  {d_core}/observed_features_correlation.qzv
        # """
        #         j_core.add_command(cmd)

        # Alpha 多样性
        d_alpha = os.path.join(dir_out, "Alpha")
        j_alpha = job.child("Alpha", d_alpha, run_type="single_run", _type="Job")
        alpha_content = ["sobs", "simpson", "shannon", "pd", "goods_coverage", "chao", "ace"]
        j_alpha.set_maxjob(len(alpha_content) + 1)

        cmd = f"""
{self.rscript} {self.bin}/rarefaction.R -d min -r \\
        -i {d_qiime2}/Denois/rename/feature.table.xls \\
        -t {d_qiime2}/Tree/rooted/tree.nwk \\
        -o {d_alpha}
{self.rscript} {self.bin}/summary_alpha.R {d_alpha}/summary.alpha_diversity.xls {d_prepare}/group.qiime.tsv {d_alpha}
"""
        j_alpha.add_command(cmd)
        for i in alpha_content:
            cmd = f"""{self.rscript} {self.bin}/smooth.R {d_alpha}/{i}.rarefaction.tsv {d_alpha}/{i}.samples.xls
{self.rscript} {self.bin}/norm2gro.R {d_alpha}/{i}.samples.xls {d_prepare}/group.qiime.tsv {d_alpha}/{i}.groups.xls
{self.generalplot} line -file {d_alpha}/{i}.samples.xls -outprefix {d_alpha}/{i}.samples
{self.generalplot} line -file {d_alpha}/{i}.groups.xls -outprefix {d_alpha}/{i}.groups
"""
            j_alpha.add_command(cmd)

        # for i in alpha_content:
        #     cmd = f"{self.env_start(self.env['qiime2'])};{self.qiime} diversity alpha --i-table {d_qiime2}/Denois/rename/feature.qza --p-metric {i} --o-alpha-diversity {d_alpha}/{i}.qza;{self.env_end()}"
        #     j_alpha.add_command(cmd)
        # cmd = f"{self.env_start(self.env['qiime2'])};{self.qiime} diversity alpha-phylogenetic --i-table {d_qiime2}/Denois/rename/feature.qza --i-phylogeny {d_tree}/rooted-tree.qza --p-metric faith_pd --o-alpha-diversity {d_alpha}/faith_pd.qza;{self.env_end()}"
        # j_alpha.add_command(cmd)

    def beta(self):
        """
        Beta diversity analysis
        """
        if (len(self.samples) >= 3):
            dir_out = os.path.join(self.out, "Beta")
            job = self.child("Beta", dir_out)

            d_prepare = os.path.join(self.out, "Prepare")
            d_qiime2 = os.path.join(self.out, "QIIME2")
            d_tree = os.path.join(d_qiime2, "Tree")
            f_tree = os.path.join(d_tree, "rooted/tree.nwk")
            f_featuretable = os.path.join(d_qiime2, "Denois/rename/feature.table.xls")

            # dist
            # TODO: 后续确认一些这里应该使用哪个表达量
            d_dist = os.path.join(dir_out, "dist")
            j_dist = job.child("dist", d_dist, _type="Job")
            cmd = f"""
{self.rscript} {self.bin}/dist.R {f_featuretable} {f_tree} {d_dist}
{self.rscript} {self.bin}/dist.R {f_featuretable} bray {d_dist}
{self.rscript} {self.bin}/dist.R {f_featuretable} jaccard {d_dist}
"""
            j_dist.add_command(cmd)

            # heatmap
            # TODO: 使用GeneralPlot内的脚本替换
            d_heatmap = os.path.join(dir_out, "heatmap")
            j_heatmap = job.child("heatmap", d_heatmap, _type="Job")
            cmd = f"""
{self.rscript} {self.bin}/heatmap.R -i {d_dist}/unweighted_unifrac.matrix.txt -o {d_heatmap}/unweighted_unifrac.heatmap -s none -c none -t "Unweighted Unifrac Heatmap"
{self.rscript} {self.bin}/heatmap.R -i {d_dist}/weighted_unifrac.matrix.txt -o {d_heatmap}/weighted_unifrac.heatmap -s none -c none -t "Weighted Unifrac Heatmap"
{self.rscript} {self.bin}/heatmap.R -i {d_dist}/bray.matrix.txt -o {d_heatmap}/bray.heatmap -s none -c none -t "Bray–Curtis Heatmap"
{self.rscript} {self.bin}/heatmap.R -i {d_dist}/jaccard.matrix.txt -o {d_heatmap}/jaccard.heatmap -s none -c none -t "Jaccard Heatmap"
"""
            j_heatmap.add_command(cmd)

            # UPGMA
            d_taxon = os.path.join(self.out, "Taxon/TaxonStat/TaxonProfile")
            d_upgma = os.path.join(dir_out, "UPGMA")
            j_upgma = job.child("upgma", d_upgma, _type="Job")
            cmd = f"""
cd {d_upgma}
for i in Phylum Class Order Family Genus Species;
do
  for j in bray jaccard weighted_unifrac unweighted_unifrac;
  do
    {self.rscript} {self.bin}/upgma_plot.R {d_dist}/$j.nwk {d_taxon}/profiling.$i.top10.xls $i $j 20
  done
done
"""
            j_upgma.add_command(cmd)

            # PCA
            # TODO: 不同的scale方法的影响
            f_group = os.path.join(d_prepare, "group.qiime.tsv")
            d_pca = os.path.join(dir_out, "PCA")
            j_pca = job.child("pca", d_pca, _type="Job")
            cmd = f"""
{self.rscript} {self.bin}/oadraw.R --scale -i {f_featuretable} --groups {f_group} -o {d_pca} -p {self.software['GeneralPlot']}
"""
            j_pca.add_command(cmd)

            # PCoA
            d_pcoa = os.path.join(dir_out, "PCoA")
            j_pcoa = job.child("pcoa", d_pcoa, _type="Job")
            cmd = f"""
for i in bray jaccard unweighted_unifrac weighted_unifrac;
do 
    {self.rscript} {self.bin}/oadraw.R --scale -t PCoA -i {f_featuretable} -d {d_dist}/$i.matrix.txt --groups {f_group} -o {d_pcoa}/$i -p {self.software['GeneralPlot']} 
done
"""
            j_pcoa.add_command(cmd)

            # NMDS
            d_nmds = os.path.join(dir_out, "NMDS")
            j_nmds = job.child("nmds", d_pcoa, _type="Job")
            cmd = f"""
for i in bray jaccard unweighted_unifrac weighted_unifrac;
do 
  {self.rscript} {self.bin}/oadraw.R --scale -t NMDS -i {f_featuretable} -d {d_dist}/$i.matrix.txt --groups {f_group} -o {d_nmds}/$i -p {self.software['GeneralPlot']} 
done
"""
            j_nmds.add_command(cmd)

    def ttest(self):
        """
        T test analysis
        """
        dir_out = os.path.join(self.out, "TTest")

        pass

    def lefse(self):
        """
        LefSe analysis
        :return:
        """
        if self.group_diff:
            dir_out = os.path.join(self.out, "LefSe")
            job = self.child("LefSe", dir_out, run_type="multi_run", _type="Job")
            job.set_maxjob(self.config["parallel"])

            d_prepare = os.path.join(self.out, "Prepare")
            d_taxon = os.path.join(self.out, "Taxon")
            for i in self.group_diff:
                compare = "-VS-".join(i)
                cmd = f"{self.python} {self.bin}/LefSe.py -i {d_taxon}/feature.taxonomy.relative/feature-table.format.xls -c {compare} -g {d_prepare}/group.qiime.tsv -o {dir_out}"
                job.add_command(cmd)

    def picrust2(self):
        """
        Function predict with PICRUSt2 analysis
        """
        dir_out = os.path.join(self.out, "PICRUSt")
        job = self.child("picrust2", dir_out, _type="Job")

        d_qiime2 = os.path.join(self.out, "QIIME2")
        f_feature = os.path.join(d_qiime2, "Denois/rename/feature.qza")
        f_seq = os.path.join(d_qiime2, "Denois/rename/rep-seqs.qza")
        job.add_command(self.env_start(self.env["qiime2"]))
        cmd = f"""set +e
{self.qiime} picrust2 full-pipeline \\
        --i-table {f_feature} \\
        --i-seq {f_seq} \\
        --p-placement-tool sepp \\
        --p-threads {self.config["threads"]} \\
        --p-max-nsti 2 \\
        --o-ko-metagenome {dir_out}/ko_metagenome.qza \\
        --o-ec-metagenome {dir_out}/ec_metagenome.qza \\
        --o-pathway-abundance {dir_out}/pathway_abundance.qza
# KO
{self.qiime} tools export \\
        --input-path {dir_out}/ko_metagenome.qza \\
        --output-path {dir_out}/ko
{self.biom} convert -i {dir_out}/ko/feature-table.biom -o {dir_out}/ko/feature-table.tsv --to-tsv
sed 1d {dir_out}/ko/feature-table.tsv | sed "1s,#OTU ID,KO,g" > {dir_out}/ko/KO.data
add_descriptions.py -i {dir_out}/ko/KO.data -m KO -o {dir_out}/ko/KO.xls

# EC
{self.qiime} tools export \\
        --input-path {dir_out}/ec_metagenome.qza \\
        --output-path {dir_out}/ec
{self.biom} convert -i {dir_out}/ec/feature-table.biom -o {dir_out}/ec/feature-table.tsv --to-tsv
sed 1d {dir_out}/ec/feature-table.tsv | sed "1s,#OTU ID,EC,g" > {dir_out}/ec/EC.data
add_descriptions.py -i {dir_out}/ec/EC.data -m EC -o {dir_out}/ec/EC.xls

# Pathway
{self.qiime} tools export \\
        --input-path {dir_out}/pathway_abundance.qza \\
        --output-path {dir_out}/pathway
{self.biom} convert -i {dir_out}/pathway/feature-table.biom -o {dir_out}/pathway/feature-table.tsv --to-tsv
sed 1d {dir_out}/pathway/feature-table.tsv | sed "1s,#OTU ID,METACYC,g" > {dir_out}/pathway/Pathway.data
add_descriptions.py -i {dir_out}/pathway/Pathway.data -m METACYC -o {dir_out}/pathway/Pathway.xls
exit 0
"""
        job.add_command(cmd)

    def package(self):
        """
        Package the result
        """
        dir_out = os.path.join(self.out, "Upload")
        job = self.child("Upload", dir_out)

        # QC
        d_i_fastqc = os.path.join(self.out, "QC/FastQC")
        d_o_fastqc = os.path.join(dir_out, "QC/FastQC")
        j_link = job.child("link", dir_out, _type="Job")
        j_link.add_command(f"mkdir -p {d_o_fastqc}")
        for sample in self.samples:
            for i in [1, 2]:
                cmd = f"ln -sf {d_i_fastqc}/{sample}_{i}_fastqc/Images/per_base_quality.png {d_o_fastqc}/{sample}_{i}.png"
                j_link.add_command(cmd)
        # ASV
        d_i_asv = os.path.join(self.out, "QIIME2/Denois/rename")
        cmd_asv = f"""
mkdir -p {dir_out}/ASV
ln -sf {d_i_asv}/dna-sequences.fasta {dir_out}/ASV/sequence.fasta
ln -sf {d_i_asv}/feature.table.xls {dir_out}/ASV/feature.table.xls
ln -sf {d_i_asv}/../stat/stats.tsv {dir_out}/ASV/denoising_stats.xls
ln -sf {d_i_asv}/../stat/denois.fill.png {dir_out}/ASV/denois.fill.png
ln -sf {d_i_asv}/../stat/denois.stack.png {dir_out}/ASV/denois.stack.png
ln -sf {d_i_asv}/../../Tree/rooted/tree.nwk {dir_out}/ASV/tree.nwk
"""
        j_link.add_command(cmd_asv)

        # taxon
        d_i_taxon = os.path.join(self.out, "Taxon")
        d_o_taxon = os.path.join(dir_out, "Taxon")
        cmd_taxon = f"""
mkdir -p {d_o_taxon}
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/all.taxonomy.stat.xls {d_o_taxon}/all.taxonomy.stat.xls
ln -sf {d_i_taxon}/TaxonStat/all.feature.tab {d_o_taxon}/feature.table.taxonomy.xls
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/taxonomy.fill.png {d_o_taxon}/taxonomy.fill.png
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/taxonomy.stack.png {d_o_taxon}/taxonomy.stack.png
mkdir -p {d_o_taxon}/profile
mkdir -p {d_o_taxon}/stack
mkdir -p {d_o_taxon}/heatmap
"""
        j_link.add_command(cmd_taxon)
        levels = ["Domain", "Phylum", "Order", "Family", "Genus", "Species"]
        for i in levels:
            cmd = f"""# {i} Stack Plot
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/profiling.{i}.xls {d_o_taxon}/profile/{i}.xls
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/profiling.{i}.top10.xls  {d_o_taxon}/stack/{i}.plot.xls
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/{i}.stack.png {d_o_taxon}/stack/{i}.stack.png
"""
            j_link.add_command(cmd)

        for i in levels[1:]:
            cmd = f"""# {i} Heatmap
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/profiling.{i}.min0.1.xls {d_o_taxon}/heatmap/{i}.plot.xls
ln -sf {d_i_taxon}/TaxonStat/TaxonProfile/{i}.heatmap.png {d_o_taxon}/heatmap/{i}.heatmap.png
"""
            j_link.add_command(cmd)

        # Diff
        if self.group_diff:
            d_i_lefse = os.path.join(self.out, "LefSe")
            d_o_lefse = os.path.join(dir_out, "Diff/LefSe")
            os.makedirs(d_o_lefse, exist_ok=True)
            for i in self.group_diff:
                compare = "-VS-".join(i)
                cmd = f"""# {compare} LefSe
ln -sf {d_i_lefse}/{compare}.Lefse.Diff.xls {d_o_lefse}/{compare}.Lefse.Diff.xls
ln -sf {d_i_lefse}/{compare}.Lefse.png {d_o_lefse}/{compare}.Lefse.png
ln -sf {d_i_lefse}/{compare}.Lefse.cladogram.png {d_o_lefse}/{compare}.Lefse.cladogram.png
"""
                j_link.add_command(cmd)

        # Alpha
        d_i_alpha = os.path.join(self.out, "Diversity/Alpha")
        d_o_alpha = os.path.join(dir_out, "Alpha")
        os.makedirs(d_o_alpha, exist_ok=True)
        alpha_content = ["sobs", "shannon", "simpson", "chao", "ace", "goods_coverage", "pd"]
        cmd = f"""
ln -sf {d_i_alpha}/summary.alpha_diversity.xls {d_o_alpha}/summary.alpha_diversity.xls
ln -sf {d_i_alpha}/summary.alpha_diversity.gro.xls {d_o_alpha}/summary.alpha_diversity.gro.xls
ln -sf {d_i_alpha}/summary.alpha_diversity.sd.xls {d_o_alpha}/summary.alpha_diversity.sd.xls
"""
        j_link.add_command(cmd)
        for i in alpha_content:
            cmd = f"""# Alpha {i}
ln -sf {d_i_alpha}/{i}.samples.png {d_o_alpha}/{i}.png
ln -sf {d_i_alpha}/{i}.groups.png {d_o_alpha}/{i}.groups.png
"""
            j_link.add_command(cmd)

        # Beta
        d_i_beta = os.path.join(self.out, "Beta")
        d_o_beta = os.path.join(dir_out, "Beta")
        if (len(self.samples) >= 3):
            # 小于三个样本无法做此类分析
            os.makedirs(d_o_beta, exist_ok=True)
            cmd = f"""# Beta
mkdir -p {d_o_beta}/heatmap
ln -sf {d_i_beta}/heatmap/*.png {d_o_beta}/heatmap
mkdir -p {d_o_beta}/UPGMA
ln -sf {d_i_beta}/UPGMA/*png {d_o_beta}/UPGMA
mkdir -p {d_o_beta}/PCA
ln -sf {d_i_beta}/PCA/*png {d_o_beta}/PCA
for i in bray jaccard unweighted_unifrac weighted_unifrac
do
  mkdir -p {d_o_beta}/PCoA/$i
  ln -sf {d_i_beta}/PCoA/$i/*png {d_o_beta}/PCoA/$i
done
for i in bray jaccard unweighted_unifrac weighted_unifrac
do
  mkdir -p {d_o_beta}/NMDS/$i
  ln -sf {d_i_beta}/NMDS/$i/*png {d_o_beta}/NMDS/$i
done
"""
            j_link.add_command(cmd)

        # PICRUSt2
        # 有些分析会没有这个结果
        d_i_picrust = os.path.join(self.out, "PICRUSt")
        d_o_picrust = os.path.join(dir_out, "PICRUSt2")
        if os.path.exists(f"{d_i_picrust}/ko/KO.xls"):
            os.makedirs(d_o_picrust, exist_ok=True)
            cmd = f"""# PICRUSt2
ln -sf {d_i_picrust}/ko/KO.xls {d_o_picrust}/KO.xls
ln -sf {d_i_picrust}/ec/EC.xls {d_o_picrust}/EC.xls
ln -sf {d_i_picrust}/pathway/Pathway.xls {d_o_picrust}/Pathway.xls
"""
            j_link.add_command(cmd)

        # Report
        j_report = job.child("report", dir_out, _type="Job")
        cmd = f"{self.perl} {self.bin}/report.pl -conf {self.f_config} -o {dir_out}"
        j_report.add_command(cmd)

        # package
        j_package = job.child('package', dir_out, _type="Job")
        cmd = f"cd {self.out}\nln -sf Upload {self.name}\nzip -r {self.name}.zip {self.name} >/dev/null 2>&1"
        j_package.add_command(cmd)

    def generate_manifest(self, samples, d_data, f_out, data_type="PE"):
        """
        Generate manifest.tsv file for QIIME2

        :param samples: The sample names to analysis
        :param d_data: The data dir
        :param f_out: The out put file name
        :param data_type: The data type(PE|SE)
        """
        if data_type == "PE":
            header = ["sample-id", "forward-absolute-filepath",
                      "reverse-absolute-filepath"]
            with open(f_out, 'w') as OUT:
                print(*header, sep="\t", file=OUT)
                for i in samples:
                    f_forward = os.path.join(d_data, f"{i}_1.fq.gz")
                    f_reverse = os.path.join(d_data, f"{i}_2.fq.gz")
                    print(*[i, f_forward, f_reverse], sep="\t", file=OUT)
        elif data_type == "SE":
            header = ["sample-id", "absolute-filepath"]
            with open(f_out, 'w') as OUT:
                print(*header, sep="\t", file=OUT)
                for i in samples:
                    f_fq = os.path.jion(d_data, f"{i}_1.fa.gz")
                    print(*[i, f_fq], sep="\t", file=OUT)
        else:
            LOG.error(f"WRONG DATA TYPE {data_type}")
            sys.exit(1)

    def generate_metadata(self, f_metadata):
        """
        Generate metadata.tsv for the qiime2 pipe

        :param f_metadata: The metadata.tsv file path
        """
        with open(f_metadata, 'w') as OUT:
            print(*["#SampleID", "Group"], sep="\t", file=OUT)
            print(*["#q2:types", "categorical"], sep="\t", file=OUT)
            for sample, group in self.sample2group.items():
                print(*[sample, group], sep="\t", file=OUT)

    def generate_qiimegroup(self, f_group):
        """
        Generate metadata.tsv for the qiime2 pipe

        :param f_group: The group.qiime.tsv file path
        """
        with open(f_group, 'w') as OUT:
            print(*["#SampleID", "Group"], sep="\t", file=OUT)
            for sample, group in self.sample2group.items():
                print(*[sample, group], sep="\t", file=OUT)

    def env_start(self, conda_env):
        """
        conda env activate command

        :param conda_env: The conda_env path
        """
        res = f"source {self.software['activate']} {conda_env}"
        return res

    def env_end(self):
        """
        conda env deactivate command
        """
        res = f"source {self.software['activate']}"
        return res
