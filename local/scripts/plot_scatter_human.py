#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.stats import pearsonr
import pdb

print('Running plot_scatter_human.py ...')

def clean_name(name):
    if '_' in name:
        name = name.split('_')[0]
        #name = name[0].upper() + name[1:]
        #name = r'\textit{%s}'%name.upper()
        name = name.upper()
    return name

sys.path.append('src')
plt.close('all')
import utils

in_ht_file = 'output/ratios_9nt_ss_locus.txt'
in_human_file = 'data/human_splice_sites.xlsx'
out_file_pattern = 'plots/psi_vs_human_scatter_%d.png'

# Load experimental psi data frame and normalize (but don't smush)
psi_df = pd.read_csv(in_ht_file, sep='\t')
psi_df.set_index('ss', inplace=True, drop=True)
cons_seq = 'CAGGUAAGU'
for col in psi_df.columns:
    psi_df.loc[:, col] = 100 * psi_df.loc[:, col] / psi_df.loc[cons_seq, col]

# Load human splice sites
tmp_df = pd.read_excel(in_human_file, header=0)

# Extract GT splice sites
human_GT_df = tmp_df.iloc[:, 0:2]
human_GT_df.columns = ['human_ct', 'seq']
human_GT_df.dropna(how='any', axis=0, inplace=True)
human_GT_df.set_index('seq', drop=True, inplace=True)
human_GT_df['GU'] = True

# Extract GC splice sites
human_GC_df = tmp_df.iloc[:, 6:8]
human_GC_df.columns = ['human_ct', 'seq']
human_GC_df.dropna(how='any', axis=0, inplace=True)
human_GC_df.set_index('seq', drop=True, inplace=True)
human_GC_df['GU'] = False

# Concatenate GT and GC splice site data frames
human_df = pd.concat([human_GT_df, human_GC_df], axis=0, join='outer')
human_df['ss'] = [ss.replace('T', 'U') for ss in human_df.index]
human_df.set_index('ss', inplace=True, drop=True)

# Merge human data with measurements
all_df = psi_df.merge(human_df, left_index=True, right_index=True, how='outer')
all_df.fillna(0, inplace=True)

# Remove cryptic splice sites
all_df = utils.remove_cryptic_splice_sites(all_df)

# Plot scatter
xcols = ('brca2_9nt', 'smn1_9nt', 'ikbkap_9nt')
ycol = 'human_ct'
xlim = [0, 150]
ylim = [-.75,1E4]
xticks = np.linspace(xlim[0], xlim[1], 4)
xticklabels = ['%d' % x for x in xticks]
fontsize = 9

for threshold in [20, 50]:
    fig = plt.figure(figsize=[8, 3.25])

    active_but_absent_df = pd.DataFrame()
    inactive_but_present_df = pd.DataFrame()
    for xcol in xcols:
        active_but_absent_df[xcol] = (all_df[xcol] > threshold) & (all_df['human_ct'] == 0)
        inactive_but_present_df[xcol] = (all_df[xcol] <= threshold) & (all_df['human_ct'] > 0)
    active_but_absent_df = active_but_absent_df[active_but_absent_df.any(axis=1)]
    active_but_absent_df.to_csv('output/active_but_absent_%d.txt'%threshold,sep='\t')
    inactive_but_present_df = inactive_but_present_df[inactive_but_present_df.any(axis=1)]
    inactive_but_present_df.to_csv('output/inactive_but_present_%d.txt' % threshold, sep='\t')

    for i, xcol in enumerate(xcols):

        # Compute histogram
        ax = fig.add_subplot(1, 3, 1+i)

        indices = np.isfinite(all_df[xcol]) & np.isfinite(all_df[ycol])
        xs = all_df[xcol][indices]
        ys = all_df[ycol][indices]

        # Compute number of active but absent sites
        num_active_and_absent = sum(active_but_absent_df[xcol])
        num_inactive_and_present = sum(inactive_but_present_df[xcol])

        # Compute correlation coefficient
        tmp_vals = ys.copy()
        tmp_vals[tmp_vals==0] = 10.0**(-0.5) # Regularization is used only for computing pearsonr
        reg_counts = np.log10(tmp_vals)
        r, pval = pearsonr(xs, reg_counts)

        # Make scatter plot
        color = 'slategrey'
        ys_jitter = ys + .67*(np.random.rand(len(ys)) - .5) #.15*np.random.randn(len(ys))
        plt.plot(xs, ys_jitter, 'o', markeredgecolor=None,
                 color=color, label='_nolegend_', alpha=.3, markersize=2)
        plt.yscale('symlog')
        plt.ylim(ylim)

        # Title
        title = 'Counts vs. %s' % (clean_name(xcol))
        plt.text(145, ylim[0]/2, 'n = %d' % num_active_and_absent, fontsize=fontsize, horizontalalignment='right', verticalalignment='center')
        #plt.text(threshold-10, 5E3, 'n=%d' % num_inactive_and_present, fontsize=fontsize, horizontalalignment='left', rotation=90)
        ax.set_title(title, fontsize=fontsize)

        # Label y axis
        ax.set_ylim(ylim)
        if i == 0:
            ax.set_ylabel('no. sites + jitter', fontsize=fontsize)
        else:
            ax.set_yticklabels([])
        ax.axhline(.5, linestyle=':', color='blue')
        ax.axvline(threshold, linestyle=':', color='tomato')

        # Label x axis
        ax.set_xlim(xlim)
        ax.set_xticks(xticks)
        ax.set_xlabel('%s (PSI)' % clean_name(xcol), fontsize=fontsize)
        ax.set_xticklabels(xticklabels, fontsize=fontsize)

    plt.tight_layout()
    out_file = out_file_pattern % threshold
    plt.savefig(out_file)
    print('Plot saved to %s' % out_file)

print('Done!')
