#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import sys

print('Running plot_histograms.py ...')

sys.path.append('src')
plt.close('all')
import utils

in_ht_file = 'output/ratios_9nt_ss_locus.txt'
out_file_format = 'plots/psi_hist_%s.pdf'

psi_df = pd.read_csv(in_ht_file, sep='\t')
psi_df.set_index('ss', inplace=True, drop=True)

# Remove cryptic splice sites
psi_df = utils.remove_cryptic_splice_sites(psi_df)

# Normalize by consensus sequence
cons_seq = 'CAGGUAAGU'
for col in psi_df.columns:
    psi_df.loc[:, col] = 100 * psi_df.loc[:, col] / psi_df.loc[cons_seq, col]
psi_df[psi_df > 100] = 100

# Add GU indicator column
psi_df['GU'] = [ss[3:5] == 'GU' for ss in psi_df.index]

for subset in ['all', 'GU', 'GC']:

    if subset == 'GU':
        tmp_df = psi_df[psi_df['GU']].copy()
        facecolor = 'steelblue'
    elif subset == 'GC':
        tmp_df = psi_df[~psi_df['GU']].copy()
        facecolor = 'lightslategrey'
    else:
        tmp_df = psi_df.copy()
        facecolor = 'midnightblue'

    fig = plt.figure(figsize=[9, 3])

    num_bins = 20
    threshold = 20
    ylim = [0, 1000]
    xlim = [0, 100]
    yticks = np.linspace(ylim[0], ylim[1], 6)
    xticks = np.linspace(xlim[0], xlim[1], 6)
    xticklabels = ['%d%%' % x for x in xticks]
    yticklabels = ['%d' % y for y in yticks]
    fontsize = 9
    cols = ['brca2_9nt', 'smn1_9nt', 'ikbkap_9nt']
    for m, col in enumerate(cols):

        # Compute histogram
        ax = fig.add_subplot(1, 3, m + 1)

        vals = tmp_df[col][np.isfinite(tmp_df[col])]
        num_sites = sum(vals > threshold)
        counts, bin_edges, xxx = \
            plt.hist(vals,
                     bins=num_bins,
                     facecolor=facecolor,
                     range=[0, 100])

        name = col.split('_')[0]
        name = name.upper() #name[0].upper() + name[1:]

        # Indicate cut-off histogram bars
        dx = bin_edges[1] - bin_edges[0]
        dy = .05 * (ylim[1] - ylim[0])
        xlb = bin_edges[0]
        xub = bin_edges[1]
        yub = ylim[1] - dy
        ylb = ylim[1] - 2 * dy
        yabove = ylim[1] + dy
        white_verts = []
        for n, x in enumerate(bin_edges[:-1]):
            if counts[n] > ylb:
                # White masking ploygon
                verts = np.array([
                    (x, ylb),
                    (x + dx, yub),
                    (x + dx, yabove),
                    (x, yabove)
                ])
                patch = patches.Polygon(verts, closed=True, facecolor='white', linewidth=1, edgecolor='white')
                ax.add_patch(patch)

                # Black line
                plt.plot((x, x + dx), (ylb, yub), linewidth=1, color='k')

        # Show threshold
        plt.axvline(threshold, linewidth=1, color='darkgray', linestyle='--')

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=fontsize)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=fontsize)
        title = '%s, %s sites\n(%d with PSI > %d%%)' % \
                (name, subset, num_sites, threshold)
        plt.title(title, fontsize=fontsize)
        plt.xlabel('PSI', fontsize=fontsize)
        plt.ylabel('splice sites', fontsize=fontsize)

    plt.tight_layout()
    out_file = out_file_format % subset
    plt.savefig(out_file)
    print('Plot saved to %s' % out_file)

print('Done!')

