#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import pandas as pd
import numpy as np

import os
import sys
import re
import pdb
import glob

print('Running analyze_efficiency.py ...')

sys.path.append('src')
import utils

# Specify input and output files
in_files = glob.glob('from_pipeline/efficiency.*.txt')
out_file_dict = {
    'flat':'output/efficiency_flat.txt',
    'num_reads':'output/efficiency_by_stage_num_reads.txt',
    'filtering':'output/efficiency_by_stage_filtering.txt'
}

# Load efficiency files
assert len(in_files)>0, 'Error! Couldnt find any input files'

# Load efficiency data into a single data frame
efficiency_df = pd.DataFrame()
for f in in_files:
    tmp_df = pd.read_csv(f,sep='\t')
    efficiency_df = efficiency_df.append(tmp_df,ignore_index=True)

# Compute useful totals.
tmp_df = efficiency_df[['LID','raw_reads']].groupby('LID').median()
print('Total number of raw_reads: %.2e'%tmp_df.sum())
print('Total number of sorted_reads: %.2e'%efficiency_df['sorted_reads'].sum())
print('Total number of parsed_reads: %.2e' % efficiency_df['parsed_reads'].sum())

# Save efficiency df
efficiency_df.to_csv(out_file_dict['flat'],sep='\t')
print('Output written to %s'%out_file_dict['flat'])

# Reformat data to num reads as one column and stage as a second column
df = pd.DataFrame()
for col in ['raw_reads','sorted_reads','parsed_reads']:
    tmp_df = pd.DataFrame()
    tmp_df['label'] = efficiency_df['group'].astype(str) +\
        ', LID ' + efficiency_df['LID'].astype(str) +\
         [(' (2x)' if 297460==lid else '') for lid in efficiency_df['LID']] 
    tmp_df['num_reads'] = efficiency_df[col]
    tmp_df['stage']=col.split('_')[0]
    df = df.append(tmp_df,ignore_index=True)

df.to_csv(out_file_dict['num_reads'],sep='\t')
print('Output written to %s'%out_file_dict['num_reads'])

# Reformat data to frac as one column and stage as a second column
df = pd.DataFrame()
for col in ['frac_sorted','frac_parsed']:
    tmp_df = pd.DataFrame()
    tmp_df['label'] = efficiency_df['group'].astype(str) +\
        ', LID ' + efficiency_df['LID'].astype(str) +\
         [(' (2x)' if 297460==lid else '') for lid in efficiency_df['LID']] 
    tmp_df['frac'] = efficiency_df[col]
    tmp_df['stage']=col.split('_')[1]
    df = df.append(tmp_df,ignore_index=True)

df.to_csv(out_file_dict['filtering'],sep='\t')
print('Output written to %s'%out_file_dict['filtering'])

print('Done!')