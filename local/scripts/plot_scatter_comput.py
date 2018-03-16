#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.stats import pearsonr
import pdb

print('Running plot_scatter_comput.py ...')

sys.path.append('src')
import utils


def clean_name(name):
    if '_' in name:
        name = name.split('_')[0]
        name = name.upper() #name[0].upper() + name[1:]
    return name


sys.path.append('src')
plt.close('all')
import utils

in_ht_file = 'output/ratios_9nt_ss_locus.txt'
out_file = 'plots/comp_vs_exp_scatter.png'

# Load experimental psi data frame and normalize (but don't smush)
psi_df = pd.read_csv(in_ht_file, sep='\t')
psi_df.set_index('ss', inplace=True, drop=True)
cons_seq = 'CAGGUAAGU'
for col in psi_df.columns:
    psi_df.loc[:, col] = 100 * psi_df.loc[:, col] / psi_df.loc[cons_seq, col]

# Remove GC splice sites
keep_indices = [ss[3:5] == 'GU' for ss in psi_df.index]
psi_df = psi_df[keep_indices]

# Remove cryptic splice sites
psi_df = utils.remove_cryptic_splice_sites(psi_df)

# Load computational predictions
comp_file = 'data/scored_splice_sites.txt'
comp_df = pd.read_csv(comp_file,
                      usecols=[0, 2, 4, 6, 8],
                      names=['ss', 'MaxEnt', 'MDD', 'MM', 'WMM'],
                      delim_whitespace=True)
comp_df['ss'] = [ss.replace('T', 'U') for ss in comp_df['ss']]

# Add in RNAhybridrid predictions
RNAhybridrid_file = 'output/rnahybrid_predictions.txt'
tmp_df = pd.read_csv(RNAhybridrid_file, sep='\t', names=['ss', 'RNAhybrid'])
comp_df = comp_df.merge(tmp_df, on='ss', how='outer')
comp_df.set_index('ss', drop=True, inplace=True)

# Remove cryptic splice sites
comp_df = utils.remove_cryptic_splice_sites(comp_df)

# Merge all into one data frame
all_df = psi_df.merge(comp_df, left_index=True, right_index=True, how='outer')

# Smush predictions
all_df = utils.smush_predictions(all_df)

xcols = ('MaxEnt', 'MDD', 'MM', 'WMM', 'RNAhybrid')
ycols = ('brca2_9nt', 'smn1_9nt', 'ikbkap_9nt')
ylim = [0, 150]
yticks = np.linspace(ylim[0], ylim[1], 4)
yticklabels = ['%d' % y for y in yticks]
fontsize = 9

fig = plt.figure(figsize=[2.5 * len(xcols), 2.5 * len(ycols)])
for i, ycol in enumerate(ycols):
    for j, xcol in enumerate(xcols):

        # Compute histogram
        panel_num = 1 + j + i * len(xcols)
        ax = fig.add_subplot(len(ycols), len(xcols), panel_num)

        indices = np.isfinite(all_df[xcol]) & np.isfinite(all_df[ycol])
        xs = all_df[xcol][indices].copy()
        ys = all_df[ycol][indices].copy()

        # Compute correlation coefficient
        r, pval = pearsonr(xs, ys)
        # mi = utils.kraskov_MI(xs,ys,6) # WHY IS THIS GIVING NEGATIVE NUMBERS!?
        # print('.')

        # Make scatter plot
        color = 'slategrey'
        plt.plot(xs, ys, 'o', markeredgecolor=None, color=color, label='_nolegend_', alpha=.5, markersize=2)

        # Title
        title = '%s vs %s, $R^2 = %.02f$' % (clean_name(ycol), clean_name(xcol), r ** 2)
        # ax.set_title(title, fontsize=fontsize, position=(0.05, 0.9),
        #              horizontalalignment='left', verticalalignment='top')
        ax.set_title(title, fontsize=fontsize)

        # Label x axis
        if i == len(ycols) - 1:
            description = 'energy\n(kcal/mol)' if xcol == 'RNAhybrid' else 'score'
            ax.set_xlabel('%s %s' % (clean_name(xcol), description), fontsize=fontsize)
        else:
            ax.set_xticklabels([])

        # Label y axis
        ax.set_ylim(ylim)
        ax.set_yticks(yticks)
        if j == 0:
            ax.set_ylabel('%s (PSI)' % clean_name(ycol), fontsize=fontsize)
            ax.set_yticklabels(yticklabels, fontsize=fontsize)
        else:
            ax.set_yticklabels([])

plt.tight_layout()
plt.savefig(out_file)
print('Plot saved to %s' % out_file)

print('Done!')

# plt.close()
