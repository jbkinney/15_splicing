#!/usr/bin/env python

from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.image as mpimg

import os
import sys
import re
import pdb
import glob
from scipy.stats import pearsonr

sys.path.append('src')
import utils

print('Running plot_validation.py ...')

in_ht_file = 'output/ratios_9nt_ss_locus.txt'
in_validation_file = 'data/library_validations.xlsx'
out_plot_file = 'plots/validation_scatter.pdf'
out_csv_file = 'output/validation.txt'

ht_df = pd.read_csv(in_ht_file, sep='\t')
ht_df.set_index('ss', inplace=True, drop=True)

# Normalize by consensus sequence
cons_seq = 'CAGGUAAGU'
for col in ht_df.columns:
    ht_df.loc[:, col] = 100 * ht_df.loc[:, col] / ht_df.loc[cons_seq, col]
ht_df[ht_df > 100] = 100

# Load validation results
lt_df = pd.read_excel(in_validation_file)
lt_df['ss'] = [ss.replace('T','U') for ss in lt_df['ss']]
lt_df.set_index('ss', inplace=True, drop=True)
lt_df.columns = [c.lower() for c in lt_df.columns]
cols_to_keep = [c for c in lt_df.columns \
                if ('validation' in c) or ('stdev' in c)]
lt_df = lt_df[cols_to_keep]

# Merge dataframes
df = ht_df.merge(lt_df, left_index=True, right_index=True, how='right')

# Remove cryptic splice sites
df = utils.remove_cryptic_splice_sites(df)

# Make plot
sns.set_style("whitegrid", {'axes.grid': False})
plt.figure(figsize=[12, 4])
loci = ['brca2', 'ikbkap', 'smn1']
for n, locus in enumerate(loci):
    plt.subplot(1, 3, n + 1)
    xs = df[locus + '_9nt']
    ys = df[locus + '_validation']
    dys = df[locus + '_stdev']
    indices = np.isfinite(xs) & np.isfinite(ys)
    r, p = pearsonr(xs[indices], ys[indices])
    plt.errorbar(xs, ys, yerr=dys, linestyle='none', marker='o')
    locus_name = locus.upper() #locus[0].upper() + locus[1:]
    plt.title('%s: $R^2$=%.2f, p=%.2e' % (locus_name, r ** 2, p))
    plt.xlabel('high throughput (PSI)')
    plt.ylabel('manual validation (PSI)')

plt.tight_layout()
plt.savefig(out_plot_file)
print('Plot saved to %s' % out_plot_file)

# Save CSV flie
df.to_csv(out_csv_file, sep='\t')
print('CSV saved to %s' % out_csv_file)

print('Done!')