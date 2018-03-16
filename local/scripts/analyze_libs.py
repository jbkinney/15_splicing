#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np

import os
import sys
import re
import pdb
import glob

print('Running analyze_libs.py ...')

sys.path.append('src')
import utils

lib_glob = 'from_pipeline/cipher.*.txt'
out_file = 'output/lib_summary.txt'

# Works: from_pipeline/cipher.brca2_9nt_ssbc_lib1.txt
pattern = re.compile('from_pipeline/cipher\.(.+_.+nt)_ssbc_(lib.+)\.txt')


# Create data frame to hold summary data
summary_df = pd.DataFrame(
    columns=['locus', 'lib',
             'num_reads',
             'num_bc',
             'num_ss',
             'frac_total_ss',
             'mean_bc_per_ss',
             'std_bc_per_ss',
             'mean_reads_per_bc',
             'std_reads_per_bc',
             'frac_good_reads']
)

# Iterate through library files

lib_files = glob.glob(lib_glob)
for i, lib_file in enumerate(lib_files):
    print('Processing %s ...'%lib_file)

    # Parse file name
    m = pattern.match(lib_file)
    assert m, 'Error! Could not parse file name %s' % lib_files
    locus = m.group(1)
    lib = m.group(2)

    # Load full bc-level information
    bc_df = pd.read_csv(lib_file, sep='\t')
    bc_df['bc_ct'] = 1

    # Keep only valid splice sites
    ss_pattern = re.compile('...G[CT]....')  # Works on CTAGCTCGA
    valid_indices = [bool(ss_pattern.match(ss)) for ss in bc_df['ss']]
    bc_df = bc_df[valid_indices]

    # Groupby splice sites
    ss_df = bc_df.groupby('ss').sum()

    if '9nt' in locus:
        max_num_ss = 2 * (4 ** 7)
    elif '11nt' in locus:
        max_num_ss = 2 * (4 ** 9)
    else:
        assert False, 'Error! Cant determine ss length for locus %s' % locus

    d = {'locus':locus, 'lib':lib}
    d['num_bc'] = len(bc_df)
    d['num_ss'] = len(ss_df)
    d['frac_total_ss'] = len(ss_df)/max_num_ss
    d['mean_bc_per_ss'] = ss_df['bc_ct'].mean()
    d['std_bc_per_ss'] = ss_df['bc_ct'].std()
    d['mean_reads_per_bc'] = bc_df['ct'].mean()
    d['std_reads_per_bc'] = bc_df['ct'].std()
    d['num_reads'] = bc_df['ct'].sum() + bc_df['otherct'].sum()
    d['frac_good_reads'] = bc_df['ct'].sum()/d['num_reads']

    summary_df.loc[i] = d

summary_df.to_csv(out_file, sep='\t')
print('Output written to %s' % out_file)

print('Done!')