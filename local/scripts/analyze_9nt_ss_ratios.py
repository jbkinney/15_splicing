#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np

import os
import sys
import re
import pdb
import glob

print('Running analyze_9nt_ss_ratios.py ...')

sys.path.append('src')
import utils

in_files_glob = 'from_pipeline/results.*_9nt_*.txt'
out_files_dict = {
    'all': 'output/ratios_9nt_ss_all.txt',
    'rep': 'output/ratios_9nt_ss_rep.txt',
    'lib': 'output/ratios_9nt_ss_lib.txt',
    'locus': 'output/ratios_9nt_ss_locus.txt',
    'psi': 'output/psi_9nt.txt',
}
reps_to_remove = ['smn1_9nt_lib1_rep1','smn1_9nt_lib3_rep3']
# Note: ratios written to file are NOT normalized to consensus

# Compute all ratios
counts_all_df = utils.get_ss_counts(in_files_glob,ss_length=9)
ratios_all_df = utils.counts_to_ratios(counts_all_df, min_tot_ct=5)
out_file = out_files_dict['all']
ratios_all_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s'%out_file)

# Save ratios with bad replicates removed
ratios_rep_df = ratios_all_df.copy()
for rep in reps_to_remove:
    if rep in ratios_rep_df.columns:
        del ratios_rep_df[rep]
out_file = out_files_dict['rep']
ratios_rep_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s'%out_file)

# Compute and save median ratios averaged over replicates in each library
tmp_df = ratios_rep_df.copy().transpose()
tmp_df['library'] = utils.splice(tmp_df.index,[0,1,2])
ratios_lib_df = tmp_df.groupby('library').median().transpose()
out_file = out_files_dict['lib']
ratios_lib_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s'%out_file)

# Compute and save median ratios averaged over all replicates for each locus
tmp_df = ratios_rep_df.copy().transpose()
tmp_df['locus'] = utils.splice(tmp_df.index,[0,1])
ratios_locus_df = tmp_df.groupby('locus').median().transpose()
out_file = out_files_dict['locus']
ratios_locus_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s'%out_file)

# Compute means, std, and rmse based on libraries
# UNLIKE ABOVE, THESE VALUES ARE NORMALIZED BY THE CONSENSUS SITE

# Define consensus sequence
cons_seq = 'CAGGUAAGU'

# Compute median PSI for each locus
psi_df = ratios_locus_df.copy()
for col in psi_df.columns:
    psi_df.loc[:, col] = \
        100 * psi_df.loc[:, col] / psi_df.loc[cons_seq, col]

# Compute median PSI for each lib
tmp_df = ratios_lib_df.copy()
for col in tmp_df.columns:
    tmp_df.loc[:, col] = \
        100 * tmp_df.loc[:, col] / tmp_df.loc[cons_seq, col]

# Compute stderr of PSI across libs at each locus
loci = ['brca2_9nt', 'ikbkap_9nt', 'smn1_9nt']
for locus in loci:
    cols = [c for c in tmp_df.columns if locus in c]
    psi_df[locus+'_stderr'] = \
        tmp_df[cols].std(axis=1) / \
        np.sqrt(np.isfinite(tmp_df[cols]).sum(axis=1))

# Order columns alphabetically
cols = list(psi_df.columns)
cols.sort()
psi_df = psi_df.loc[:, cols]

# Save PSI file
out_file = out_files_dict['psi']
psi_df.to_csv(out_file, sep='\t', na_rep='NaN')
print('Saving output to %s'%out_file)

print('Done!')