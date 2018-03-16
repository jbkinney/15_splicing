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

print('Running analyze_11nt_ss_ratios.py ...')

sys.path.append('src')
import utils

in_files_glob = 'from_pipeline/results.*_11nt_*.txt'
out_files_dict = {
    'rep':'output/ratios_11nt_ss_rep.txt',
    'lib':'output/ratios_11nt_ss_lib.txt',
    'locus':'output/ratios_11nt_ss_locus.txt',
    'rep_marg':'output/ratios_11nt_ss_rep_marg.txt',
    'lib_marg':'output/ratios_11nt_ss_lib_marg.txt',
    'locus_marg':'output/ratios_11nt_ss_locus_marg.txt'
}
# Note: ratios written to file are NOT normalized to consensus

# Compute all ratios
counts_rep_df = utils.get_ss_counts(in_files_glob,ss_length=11)
ratios_rep_df = utils.counts_to_ratios(counts_rep_df, min_tot_ct=5)
out_file = out_files_dict['rep']
ratios_rep_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s' % out_file)

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

# Marginalize counts
counts_rep_df['ss_9nt'] = [ss[:-2] for ss in counts_rep_df.index]
counts_rep_marg_df = counts_rep_df.groupby('ss_9nt').sum()
counts_rep_marg_df.index.name = 'ss'

# Compute ratios for marginalized counts
ratios_rep_marg_df = utils.counts_to_ratios(counts_rep_marg_df, min_tot_ct=5)
out_file = out_files_dict['rep_marg']
ratios_rep_marg_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s'%out_file)

# Compute and save median marg ratios averaged over reps in each library
tmp_marg_df = ratios_rep_marg_df.copy().transpose()
tmp_marg_df['library'] = utils.splice(tmp_marg_df.index,[0,1,2])
ratios_lib_marg_df = tmp_marg_df.groupby('library').median().transpose()
out_file = out_files_dict['lib_marg']
ratios_lib_marg_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s'%out_file)

# Compute and save median marg ratios averaged over all reps at each locus
tmp_marg_df = ratios_rep_marg_df.copy().transpose()
tmp_marg_df['locus'] = utils.splice(tmp_marg_df.index,[0,1])
ratios_locus_marg_df = tmp_marg_df.groupby('locus').median().transpose()
out_file = out_files_dict['locus_marg']
ratios_locus_marg_df.to_csv(out_file,sep='\t',na_rep='NaN')
print('Saving output to %s'%out_file)

print('Done!')