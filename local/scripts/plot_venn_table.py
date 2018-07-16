#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.cm as cm

import sys
sys.path.append('src')
import utils
from helper_functions import gelx, gely

in_file = 'output/ratios_9nt_ss_locus.txt'
out_file = 'plots/venn_table.pdf'


print('Running plot_venn_table.py ...')

# Load ratios
ratios_df = pd.read_csv(in_file, sep='\t')
ratios_df.set_index('ss', inplace=True, drop=True)

# Normalize by consensus sequence
psi_df = ratios_df.copy()
cons_seq = 'CAGGUAAGU'
for col in ratios_df.columns:
    psi_df.loc[:, col] = 100 * ratios_df.loc[:, col] / ratios_df.loc[cons_seq, col]
psi_df[psi_df > 100] = 100 - 1E-6

# Remove cryptic splice sites
psi_df = utils.remove_cryptic_splice_sites(psi_df)



# Set thresholds
threshold_dict = {
    'lo':(0, 20),
    'med':(20, 80),
    'hi':(80, 100)
}
categories = ['lo','med','hi']

# Create table of counts
counts_df = pd.DataFrame()
for v_col in psi_df.columns:
    for h_col in psi_df.columns:
        for v_cat in categories:
            for h_cat in categories:
                v_name = '%s_%s'%(v_col, v_cat)
                h_name = '%s_%s'%(h_col, h_cat)
                v_lb, v_ub = threshold_dict[v_cat]
                h_lb, h_ub = threshold_dict[h_cat]
                counts_df.loc[v_name,h_name] = sum( 
                    (psi_df[v_col] >= v_lb) & 
                    (psi_df[v_col] < v_ub) & 
                    (psi_df[h_col] >= h_lb) & 
                    (psi_df[h_col] < h_ub)
                )        
counts_df = counts_df.astype(int)



# Create annotation df
annotation_df = pd.DataFrame()
annotation_df['locus']=['']*9
annotation_df['category']=['']*9
row = 0
for col in psi_df.columns:
    for cat in categories:
        locus = col.split('_')[0].upper()
        annotation_df.iloc[row,:] = [locus,cat]
        row += 1
annotation_df



# Draw heatmap
fig = plt.figure(figsize=(6, 5))
left = .15
bottom = .05
width = .8
height = .8
fontsize = 9
cmap = cm.get_cmap('YlGn')
annotation_spacing = .5

ax = fig.add_axes([left, bottom, width, height],aspect='equal')
sns.heatmap(counts_df,
            ax=ax, vmax=1000, annot=True, fmt='d',
            annot_kws={"size": 7}, 
            cmap=cmap, cbar_kws={"pad": .05, "aspect": 50})
gelx(ax, annotation_df,
     fontsize=fontsize,
     annotation_spacing=annotation_spacing)
gely(ax, annotation_df,
     annotation_spacing=annotation_spacing,
     fontsize=fontsize)

for x in [2.99,5.99]:
    ax.axhline(x,linewidth=1,color='k')
    ax.axvline(x,linewidth=1,color='k')

ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
    
# Change fontsize of colorbar tick labels
cax = plt.gcf().axes[-1]
cax.tick_params(labelsize=fontsize)
fig.savefig(out_file)
print('Done! Output saved to %s'%out_file)
