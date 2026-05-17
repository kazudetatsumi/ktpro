#!/usr/bin/env python
"""
This script pre-processes a text file for latent Dirichlet Allocation,
forming pands data frame and saved as a csv file.
The original ipynb script written by Dr. Yoshifumi Amamoto.
Transferred it to this file by Kazuyoshi TATSUMI 2026/05/17.
"""
import pandas as pd

index_list = ['Reference Type', 'Year', 'Title', 'Author', 'Journal',
              'Abstract', 'Times Cited', 'ISSN']

all_series = []

with open('Neutron.txt', 'r') as f:
    df_inter = None
    for line in f:
        line = line.rstrip()
        if 'Reference Type:' in line:
            if df_inter is not None:
                all_series.append(df_inter)
            df_inter = pd.Series(index=index_list, dtype=object)
            # 【ここが重要】すべてのSeriesの名前を 0 に固定します
            # これにより転置後のインデックスがすべて 0 になり、元のバグ挙動を再現します
            df_inter.name = 0
            df_inter['Reference Type'] = line[18:]
            continue
        if df_inter is not None:
            if line[0:5] == 'Title':
                df_inter['Title'] = line[7:]
            elif line[0:7] == 'Author:':
                df_inter['Author'] = line[7:]
            elif line[0:4] == 'Year':
                df_inter['Year'] = line[6:]
            elif line[0:7] == 'Journal':
                df_inter['Journal'] = line[9:]
            elif line[0:8] == 'Abstract':
                df_inter['Abstract'] = line[10:]
            elif line[0:11] == 'Times Cited':
                df_inter['Times Cited'] = line[13:]
            elif line[0:4] == 'ISSN':
                df_inter['ISSN'] = line[6:]
if df_inter is not None:
    all_series.append(df_inter)
# axis=1 で結合すると、列名がすべて 0 の DataFrame ができます
df = pd.concat(all_series, axis=1)
# df.T をすると、行名（インデックス）がすべて 0 になります
df.T.to_csv('Neutron.csv')
