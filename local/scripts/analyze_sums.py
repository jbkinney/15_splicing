#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np

import os
import sys
import re
import pdb
import glob

print('Running analyze_sums.py ...')

# Import my stuff
sys.path.append('src')
import utils

# Define input and output files
in_9nt_bc_lib_glob = 'from_pipeline/results.*_9nt_*.txt'
in_jct_glob = 'from_pipeline/counts.*jct*.txt'
out_bc_sums_file = 'output/sums_bc.txt'
out_lib_sums_file = 'output/sums_lib.txt'
out_jct_sums_file = 'output/sums_jct.txt'

#
# Compute the total number of reads in each experiment and each library
#
print('Counting num reads in all bc and lib samples...')

# Load read counts from all experiments (except junction experiments)
counts_9nt_df = utils.get_ss_counts(in_9nt_bc_lib_glob,ss_length=9)
all_sums = counts_9nt_df.sum(axis=0)

# Sum counts
experiments = utils.splice(all_sums.index,elements=[0,1,2,3],unique=True)
cols = utils.splice(all_sums.index,elements=[4,5],unique=True)
all_sums_df = pd.DataFrame(index=experiments,columns=cols)
all_sums_df.index.name='experiment'
for experiment in experiments:
    for col in cols:
        sample = '%s_%s'%(experiment,col)
        all_sums_df.loc[experiment,col] = int(all_sums[sample])
    all_sums_df.loc[experiment,'library'] = \
        '_'.join(experiment.split('_')[:3])
        
# Get sums just for bc samples
bc_sums_df = all_sums_df[['tot_ct','ex_ct']]
bc_sums_df.to_csv(out_bc_sums_file,sep='\t')
print('Sums written to %s'%out_bc_sums_file)

# Get sums just for libraries
lib_sums_df = all_sums_df.groupby('library').first()[['lib_ct','mis_ct']]
lib_sums_df.to_csv(out_lib_sums_file,sep='\t')
print('Sums written to %s'%out_lib_sums_file)

#
# Compute the total number of reads in each junction sample
#
print('Counting num reads in all jct samples...')

file_names = glob.glob(in_jct_glob)
jct_sums_df = pd.DataFrame()
jct_sums_df.index.name='experiment'
for file_name in file_names:

    # Parse file_name
    m = re.match(utils.COUNT_FILE_PATTERN,file_name)
    d = m.groupdict()
    sample = d['sample']
    if 'brca2' in sample:
        sample += '_rep1'
    experiment = utils.splice(sample,elements=[0,1,3,4])

    # Load data
    tmp_df = pd.read_csv(file_name,delim_whitespace=True)
    
    # Marginalize over barcodes; keep only ex and tot cts. 
    jct_sums_df.loc[experiment,'jct_ct'] = tmp_df['ct'].sum()
    
jct_sums_df = jct_sums_df.astype(int)
jct_sums_df.to_csv(out_jct_sums_file,sep='\t')
print('Sums written to %s'%out_jct_sums_file)

print('Done!')