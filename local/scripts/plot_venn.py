#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib_venn import venn3

print('Running plot_venn.py ...')

sys.path.append('src')
import utils

in_file = 'output/ratios_9nt_ss_locus.txt'
out_file = 'plots/venn.pdf'

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

a_col = 'brca2_9nt'
b_col = 'smn1_9nt'
c_col = 'ikbkap_9nt'

fontsize = 9

# Format labels
set_labels = (a_col, b_col, c_col)
set_labels = [x.split('_')[0] for x in set_labels]
set_labels = [x[0].upper() + x[1:] for x in set_labels]

threshold_list = [
    (0, 20),
    (20, 80),
    (80, 100)
]

fig = plt.figure(figsize=[9, 3])
for n, (low, hi) in enumerate(threshold_list):
    ax = fig.add_subplot(1, 3, n + 1)
    venn_dict = {}
    for a in [0, 1]:
        for b in [0, 1]:
            for c in [0, 1]:
                key = '%d%d%d' % (a, b, c)
                indices = \
                    ((low <= psi_df[a_col]) & (hi > psi_df[a_col]) == bool(a)) & \
                    ((low <= psi_df[b_col]) & (hi > psi_df[b_col]) == bool(b)) & \
                    ((low <= psi_df[c_col]) & (hi > psi_df[c_col]) == bool(c))
                venn_dict[key] = indices.sum()

    # Draw Venn diagram
    venn_obj = venn3(subsets=venn_dict, set_labels=set_labels, ax=ax)
    plt.title('%d%% $\leq$ PSI < %d%%' % (low, hi), fontsize=fontsize)

    # Set label size
    for t in venn_obj.set_labels:
        if t:
            t.set_fontsize(0)
            t.set_weight('bold')
            t.set_alpha(0.0)
    for t in venn_obj.subset_labels:
        if t:
            t.set_fontsize(9)
            t.set_alpha(0.0)

plt.savefig(out_file)
print('Saving output to file %s' % out_file)
print('Done!')
