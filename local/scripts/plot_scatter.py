#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.stats import pearsonr

sys.path.append('src')
import utils

print('Running plot_scatter.py ...')

def clean_name(name):
    name = name.split('_')[0]
    name = name.upper() #name[0].upper() + name[1:]
    return name

in_ht_file = 'output/ratios_9nt_ss_locus.txt'
out_file_format = 'plots/psi_scatter_%s.png'

psi_df = pd.read_csv(in_ht_file, sep='\t')
psi_df.set_index('ss', inplace=True, drop=True)

# Remove cryptic splice sites
psi_df = utils.remove_cryptic_splice_sites(psi_df)

# Normalize by consensus sequence
cons_seq = 'CAGGUAAGU'
for col in psi_df.columns:
    psi_df.loc[:, col] = 100 * psi_df.loc[:, col] / psi_df.loc[cons_seq, col]

# Add GU indicator column
psi_df['GU'] = [ss[3:5] == 'GU' for ss in psi_df.index]

# Remove rows corresponding to cryptic sites
cryptic_site_indices = np.array([ss[1:3] == 'GU' for ss in psi_df.index])
psi_df = psi_df[~cryptic_site_indices]

special_sites = [
    ('consensus', 'CAGGUAAGU', {'color': 'deepskyblue', 'marker': 'o'}),
    ('BRCA2 WT', 'CAGGCAAGU', {'color': 'coral', 'marker': 'o'}),
    ('SMN1 WT', 'GGAGUAAGU',{'color': 'orchid', 'marker': 'o'}),
    ('IKBKAP WT', 'CAAGUAAGU', {'color': 'yellowgreen', 'marker': 'o'}),
    ('IVS17+5G>A', 'CAGGCAAAU', {'color': 'coral', 'marker': '^'}),
    ('IVS17-1G>C', 'CACGCAAGU', {'color': 'coral', 'marker': 'v'}),
    ('BRCA2 M1', 'CAGGCCAGU', {'color': 'coral', 'marker': '<'}),
    ('BRCA2 M2', 'CAGGCACGU', {'color': 'coral', 'marker': '>'}),
    ('SMN1 M1', 'GGAGUCAGU', {'color': 'orchid', 'marker': '^'}),
    ('SMN1 M2', 'GGAGUAAAU', {'color': 'orchid', 'marker': 'v'}),
    ('IVS20+6T>C', 'CAAGUAAGC', {'color': 'yellowgreen', 'marker': '^'}),
]


for subset in ['GY', 'allmut', 'all', 'GU', 'GC']:

    if subset == 'GU':
        tmp_df = psi_df[psi_df['GU']].copy()
        color = 'steelblue'
    elif subset == 'GC':
        tmp_df = psi_df[~psi_df['GU']].copy()
        color = 'teal'
    else:
        tmp_df = psi_df.copy()
        color = 'lightslategrey'

    lim = [0, 150]
    ticks = np.linspace(lim[0], lim[1], 4)
    ticklabels = ['%d' % x for x in ticks]
    fontsize = 9
    pairs = [
        ('brca2_9nt', 'smn1_9nt'),
        ('brca2_9nt', 'ikbkap_9nt'),
        ('smn1_9nt','ikbkap_9nt')]

    num_panels = 4 if 'all' in subset else 3
    fig = plt.figure(figsize=[3*num_panels, 3.25])
    if 'all' in subset:
        pairs.append(('', ''))

    for m, pair in enumerate(pairs):

        # Compute histogram
        ax = fig.add_subplot(1, num_panels, m + 1, aspect='equal')

        if m < 3:
            colx = pair[0]
            coly = pair[1]
            indices = np.isfinite(tmp_df[colx]) & np.isfinite(tmp_df[coly])
            xs = tmp_df[colx][indices]
            ys = tmp_df[coly][indices]

            # Compute correlation coefficient
            r, pval = pearsonr(xs,ys)

            # Make scatter plot
            plt.plot(xs, ys, 'o', markeredgecolor=None, color=color,
                     label='_nolegend_', alpha=.5, markersize=2)
            title = '%s vs %s\n(%s sites, $R^2 = %.02f)$' % (clean_name(coly), clean_name(colx), subset, r ** 2)

        # Plot special sites
        if 'all' in subset:
            if subset == 'all':
                sites = special_sites[:4]
            elif subset == 'allmut':
                sites = special_sites
            else:
                assert False, 'Error: dont recognize subset %s'%subset

            for label, ss, kwargs in sites:
                label = '%s (%s)' % (label, ss)
                if ss in xs.index and ss in ys.index:
                    if m < 3:
                        x = xs[ss]
                        y = ys[ss]
                    elif m == 3:
                        x = -100
                        y = -100
                    else:
                        assert False, 'Error! Invalid m=%d' % m

                    plt.plot(x, y, linestyle='none', markersize=7, label=label,
                             markeredgecolor='k', markeredgewidth=.5,
                             **kwargs)
                if m == 3:
                    plt.legend(loc='center')
                    plt.axis('off')
                    title = ''
                    # if subset == 'all':
                    #     plt.legend(bbox_to_anchor=(1.15, .75), ncol=1)
                    # elif subset == 'allmut':
                    #     plt.legend(bbox_to_anchor=(1.15, .05), ncol=1)

        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticklabels, fontsize=fontsize)
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels, fontsize=fontsize)
        plt.title(title, fontsize=fontsize)
        plt.xlabel('%s (PSI)'%clean_name(colx), fontsize=fontsize)
        plt.ylabel('%s (PSI)'%clean_name(coly), fontsize=fontsize)

    plt.tight_layout()
    out_file = out_file_format % subset
    plt.savefig(out_file)
    print('Plot saved to %s' % out_file)

print('Done!')

