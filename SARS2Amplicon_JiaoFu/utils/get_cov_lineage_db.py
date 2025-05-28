#!/usr/bin/env python
from pathlib import Path
import pandas as pd

cov_lineage_table = f"{Path(__file__).parents[1]}/etc/Cov_Lineage.txt"
cov_lineage_html = "https://cov-lineages.org/lineage_list.html"
res = pd.read_html(cov_lineage_html)
df_covlineage = res[0]
df_covlineage.fillna("-", inplace=True)
df_covlineage.to_csv(cov_lineage_table, sep="\t", index=False, encoding="utf-8")
