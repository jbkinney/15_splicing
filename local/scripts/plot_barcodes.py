#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import os
import sys
import re
import pdb
import glob

print('Running plot_barcodes.py ...')

sys.path.append('src')
plt.close('all')
import utils

# Set file names
data_set = 'brca2_9nt_lib1_rep2'
in_data_file = 'from_pipeline/results.%s.txt'%data_set
out_file = 'plots/barcode.%s.pdf'%data_set
cons_ss = 'CAGGTAAGT'

# Load data frame
print('Loading %s...'%in_data_file)
df = pd.read_csv(in_data_file, sep='\t')

# Get ratio for consensus splice site
cons_df = df[df['ss'] == cons_ss]
cons_ratio = cons_df['ex_ct'].sum() / cons_df['tot_ct'].sum()

# Filter barcodes to only include those with tot_ct > 10
df = df[df['tot_ct'] > 10]
df['ratio'] = df['ex_ct']/df['tot_ct']

# Compute median ratios
median_ratios = df[['ss','ratio']].groupby('ss').median()
median_ratios.sort_values('ratio',ascending=False, inplace=True)
median_ratios['num_bc'] = df.groupby('ss').size()
median_ratios['psi'] = 100*median_ratios['ratio'] / cons_ratio
median_ratios.head()

# Filter splice sites to include only those with
# median PSI > 20 and num_bc > 10
indices = (median_ratios['psi'] > 20) & (median_ratios['num_bc'] > 10)
print('Number of ss passing all thresholds: %d'%sum(indices))

# Keep at most 20 of these splice sites
all_ss_to_keep = median_ratios.index[indices]
ss_to_keep = np.random.choice(all_ss_to_keep, 
                              min(20,len(all_ss_to_keep)), 
                              replace=False)
tmp_df = df[df['ss'].isin(ss_to_keep)]
tmp_df['psi'] = 100 * tmp_df['ratio'] / cons_ratio

# Make box plot
sns.set_style('white')
sns.boxplot(x='ss', y='psi', data=tmp_df, color='white')
sns.stripplot(x='ss', y='psi', data=tmp_df, alpha=.5)
plt.ylim([0,200])
plt.xticks([])
sns.despine(offset=10, trim=True);
plt.title('PSI by barcode for %s'%data_set)

print('Writing boxplot to %s...'%out_file)
plt.savefig(out_file)

print('Done!')

